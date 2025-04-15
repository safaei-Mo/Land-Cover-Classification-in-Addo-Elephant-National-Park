Map.centerObject(Addo_Boundary);
var SRTM = ee.Image("CGIAR/SRTM90_V4");
var L_8_ImageCollectionSR = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2");


var Landsat_8_Collection = L_8_ImageCollectionSR
                  .filterBounds(Addo_Boundary)
                  .filterDate('2014-01-01','2014-12-31')
                  .filter(ee.Filter.lt("CLOUD_COVER",10))
                  
                  
var Landsat_8_Bands = Landsat_8_Collection
                  .select(["SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7"])
                  .median()
                  .clip(Addo_Boundary);

var temporalCollection = function(collection, start, count, interval, units){
  var sequences = ee.List.sequence(0,ee.Number(count).subtract(1.0));
  var originDate = ee.Date(start);
  
  return ee.ImageCollection(sequences.map(function(i){
    
    var startDate = originDate.advance(ee.Number(interval).multiply(i),units);
    var endDate = originDate.advance(ee.Number(interval).multiply(ee.Number(i).add(1)),units);
    
    return collection.filterDate(startDate, endDate)
    .median()
    .select(["SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","ST_B10"])
    .clip(Addo_Boundary)
    .set('systam:time_start',startDate.millis())
    .set('system:time_end', endDate.millis());
  }));
};

var NDI_Function_1 = function(img){
  return img.normalizedDifference(['SR_B2','SR_B1'])};

var NDI_Function_2 = function(img){
  return img.normalizedDifference(['SR_B3','SR_B2'])};
  
var NDI_Function_3 = function(img){
  return img.normalizedDifference(['SR_B4','SR_B3'])};
  
var NDI_Function_4 = function(img){
  return img.normalizedDifference(['SR_B5','SR_B4'])};
  
var NDI_Function_5 = function(img){
  return img.normalizedDifference(['SR_B6','SR_B5'])};
  
var NDI_Function_6 = function(img){
  return img.normalizedDifference(['SR_B7','SR_B6'])};

var Landsat_8_TemporalCollection = temporalCollection(L_8_ImageCollectionSR, '2014-01-01', 6, 2, 'month')


var NDI_1 = Landsat_8_TemporalCollection.map(NDI_Function_1).toBands()
var NDI_2 = Landsat_8_TemporalCollection.map(NDI_Function_2).toBands()
var NDI_3 = Landsat_8_TemporalCollection.map(NDI_Function_3).toBands()
var NDI_4 = Landsat_8_TemporalCollection.map(NDI_Function_4).toBands()
var NDI_5 = Landsat_8_TemporalCollection.map(NDI_Function_5).toBands()
var NDI_6 = Landsat_8_TemporalCollection.map(NDI_Function_6).toBands()
var NDI_Stacked = ee.Image.cat([NDI_1,NDI_2,NDI_3,NDI_4,NDI_5,NDI_6]);


var NDVI_Function = function(img){
  return img.normalizedDifference(['SR_B5','SR_B4'])};

var NDVI_6 = Landsat_8_TemporalCollection.map(NDVI_Function).toBands()



var computeGLCM = function(img){
  var glcm_bands = [
       'corr', 'var',
     'sent','ent','savg','diss',
    'inertia','prom'];
    
  var bandNames = img.bandNames().getInfo();
  
  var selected_glcm_bands = ee.List(bandNames.map(function(i){
    return glcm_bands.map(function(j){
      return i + '_' + j;
    });
  })).flatten();
  
  return img.glcmTexture({
    size: 4
  }).select(selected_glcm_bands);
};

var glcm = computeGLCM(PCA_Landsat_8.toInt32());
                  

var NDWI_Function = function(img){
  return img.normalizedDifference(['SR_B3','SR_B5'])};

var NDWI = Landsat_8_Collection
                  .map(NDWI_Function)
                  .median()
                  .clip(Addo_Boundary);




var SR_B6 = Landsat_8_Bands.select('SR_B6');    
var SR_B4 = Landsat_8_Bands.select('SR_B4');    
var SR_B5 = Landsat_8_Bands.select('SR_B5');    
var SR_B2 = Landsat_8_Bands.select('SR_B2');  

var B6_add_B4 = SR_B6.add(SR_B4)
var B5_add_B2 =  SR_B5.add(SR_B2)

var numerator = B6_add_B4.subtract(B5_add_B2);
var denominator = B6_add_B4.add(B5_add_B2);

var BSI = numerator.divide(denominator);
var BSI_2= BSI.multiply(100).add(100);


var nbi_Function = function(img){
  return img.expression(
    '(swir * red) / nir',
    {
        red: img.select('SR_B4'),    
        nir: img.select('SR_B5'),    
        swir: img.select('SR_B6')   
    })};

var NBI = Landsat_8_Collection
                  .map(nbi_Function)
                  .median()
                  .clip(Addo_Boundary);

var SRTM = SRTM.resample('bilinear').reproject({
  crs: SRTM.projection(),
  scale: 30
});
var elevation = SRTM.select('elevation').clip(Addo_Boundary);
var slope = ee.Terrain.slope(elevation);
var aspect = ee.Terrain.aspect(elevation);
var terrain = ee.Image.cat([elevation, slope, aspect]);



var Independetn_Variables = ee.Image.cat([NDI_Stacked,NDVI_6, glcm,NDWI,BSI_2,NBI,terrain,SMA]);

print("Independent_Variables",Independetn_Variables);


var Training_values = Independetn_Variables.sampleRegions({
collection:Training_Samples_2014,
properties:['class'],
scale:30,
geometries: true
});



var classifier = ee.Classifier.smileRandomForest({numberOfTrees:1000, variablesPerSplit:9}).train({
features:Training_values,
classProperty:'class',
inputProperties:Independetn_Variables.bandNames(), 
});

var classification_Map = Independetn_Variables.classify(classifier);



Map.addLayer(classification_Map,{min:1,max:6, palette:["blue","red","brown","green","orange","yellow"]},"classified",false);



var Test_values = Independetn_Variables.sampleRegions({
collection:Test_Samples_2014,
properties:['class'],
scale:30,
geometries: true
});


var test = Test_values.classify(classifier);

var error_matrix = test.errorMatrix('class', 'classification');
print("OA", error_matrix.accuracy());
print('F1 Score:', error_matrix.fscore(1.0));



var importance = ee.Dictionary(classifier.explain().get('importance'));
var sum = importance.values().reduce(ee.Reducer.sum());
var relativeImportance = importance.map(function(key, val) {
  return ee.Number(val).multiply(100).divide(sum);
});

print("Variable Importance", relativeImportance);
