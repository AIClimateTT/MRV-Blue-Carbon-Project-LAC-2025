/*******************************************************
METHOD 1 — PSEUDO TRAINING EXPORT (GEE) — CHM FOOTPRINT REGION
Uses CHM coverage as the sampling/export region (no AOI needed).

Exports pseudo_train.csv with Lon/Lat in EPSG:4326:
- CHM metrics per GRID cell (mean, p95, stdDev, cover_2m_mean, cover_5m_mean, valid_frac)
- S1 predictors (VV, VH, VH_VV_ratio, VV_minus_VH)
- S2 predictors (B2,B3,B4,B8,B11,B12, NDVI,EVI,NDWI,NDRE)
- tile_id (meter-based)
- Longitude, Latitude
*******************************************************/

var CHM_ASSET
var WATER_RASTER_ASSET

var START = '2022-07-02';
var END   = '2023-09-14';

var DRIVE_FOLDER = 'BlueCarbonExports';
var EXPORT_NAME  = 'Colombia_pseudo_train_grid20m_2022_2023';

var GRID_SIZE = 20;        // meters
var MAX_SAMPLES = 200000;

var VALID_FRAC_MIN = 0.80;
var CHM_P95_MIN = 0.5;
var CHM_P95_MAX = 80;

// Optional thinning if export is huge (1.0 keeps all)
var THIN_FRAC = 1.0;

// ESRI:103599 MAGNA-SIRGAS_CMT12 WKT
var WKT_103599 =
'PROJCS["MAGNA-SIRGAS_CMT12",GEOGCS["MAGNA-SIRGAS",DATUM["Marco_Geocentrico_Nacional_de_Referencia",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6686"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4686"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",4],PARAMETER["central_meridian",-73],PARAMETER["scale_factor",0.9992],PARAMETER["false_easting",5000000],PARAMETER["false_northing",2000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["ESRI","103599"]]';

var REF_PROJ = ee.Projection(WKT_103599);
var CHM_SCALE = 0.3;

// Helpers
function forceDefaultProj(img, scaleMeters) {
  return ee.Image(img).setDefaultProjection(REF_PROJ, null, scaleMeters);
}

function aggToGrid(img, reducer, bandName) {
  img = forceDefaultProj(img, CHM_SCALE);
  return img
    .reduceResolution({
      reducer: reducer,
      maxPixels: 16384
    })
    .reproject({
      crs: REF_PROJ,
      scale: GRID_SIZE
    })
    .rename(bandName);
}

/**** =========================
 *  1) LOAD CHM AND DEFINE REGION = CHM FOOTPRINT
 *  ========================= ****/

var chm = forceDefaultProj(ee.Image(CHM_ASSET), CHM_SCALE).rename('CHM');
var REGION = chm.geometry();   // ✅ use CHM footprint as region

var water = forceDefaultProj(ee.Image(WATER_RASTER_ASSET), CHM_SCALE).rename('water'); // 0 land, 1 water

/**** =========================
 *  2) CLEAN CHM
 *  ========================= ****/

var chmClean = chm
  .updateMask(water.eq(0))
  .updateMask(chm.gte(0))
  .updateMask(chm.lte(CHM_P95_MAX));

var cover2_native = forceDefaultProj(chmClean.gt(2), CHM_SCALE);
var cover5_native = forceDefaultProj(chmClean.gt(5), CHM_SCALE);
var valid_native  = forceDefaultProj(chmClean.mask(), CHM_SCALE);

/**** =========================
 *  3) CHM METRICS ON GRID
 *  ========================= ****/

var chm_mean   = aggToGrid(chmClean, ee.Reducer.mean(), 'CHM_mean');
var chm_p95    = aggToGrid(chmClean, ee.Reducer.percentile([95]), 'CHM_p95');
var chm_std    = aggToGrid(chmClean, ee.Reducer.stdDev(), 'CHM_stdDev');
var cover2     = aggToGrid(cover2_native, ee.Reducer.mean(), 'cover_2m_mean');
var cover5     = aggToGrid(cover5_native, ee.Reducer.mean(), 'cover_5m_mean');
var valid_frac = aggToGrid(valid_native, ee.Reducer.mean(), 'valid_frac');

var chmStack = ee.Image.cat([chm_mean, chm_p95, chm_std, cover2, cover5, valid_frac]);

/**** =========================
 *  4) SENTINEL-2 COMPOSITE
 *  ========================= ****/

function addS2Indices(img) {
  var scaled = img.divide(10000);
  var b2 = scaled.select('B2');
  var b3 = scaled.select('B3');
  var b4 = scaled.select('B4');
  var b5 = scaled.select('B5');
  var b8 = scaled.select('B8');

  var ndvi = b8.subtract(b4).divide(b8.add(b4)).rename('NDVI');
  var ndwi = b3.subtract(b8).divide(b3.add(b8)).rename('NDWI');
  var evi  = b8.subtract(b4).multiply(2.5)
    .divide(b8.add(b4.multiply(6)).subtract(b2.multiply(7.5)).add(1))
    .rename('EVI');
  var ndre = b8.subtract(b5).divide(b8.add(b5)).rename('NDRE');

  return scaled.addBands([ndvi, ndwi, evi, ndre]);
}

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate(START, END)
  .filterBounds(REGION)
  .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 80))
  .map(addS2Indices)
  .median()
  .select(['B2','B3','B4','B8','B11','B12','NDVI','EVI','NDWI','NDRE']);

/**** =========================
 *  5) SENTINEL-1 COMPOSITE
 *  ========================= ****/

function addS1Derived(img) {
  var vv = img.select('VV');
  var vh = img.select('VH');
  return img.addBands([
    vh.divide(vv).rename('VH_VV_ratio'),
    vv.subtract(vh).rename('VV_minus_VH')
  ]);
}

var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterDate(START, END)
  .filterBounds(REGION)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .map(addS1Derived)
  .median()
  .select(['VV','VH','VH_VV_ratio','VV_minus_VH']);

/**** =========================
 *  6) AGGREGATE SATELLITES TO GRID
 *  ========================= ****/

s1 = s1.setDefaultProjection(REF_PROJ, null, 10);
s2 = s2.setDefaultProjection(REF_PROJ, null, 10);

function satToGrid(img) {
  return ee.Image(img)
    .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 1024})
    .reproject({crs: REF_PROJ, scale: GRID_SIZE});
}

var s1g = satToGrid(s1);
var s2g = satToGrid(s2);

/**** =========================
 *  7) EXPORT IMAGE + SAMPLE
 *  ========================= ****/

var exportImg = ee.Image.cat([chmStack, s1g, s2g]).clip(REGION);

var samples = exportImg.sample({
  region: REGION,
  scale: GRID_SIZE,
  projection: REF_PROJ.atScale(GRID_SIZE),
  geometries: true,
  numPixels: MAX_SAMPLES,
  tileScale: 16
});

// Add lon/lat + tile_id
samples = samples.map(function(f) {
  var geom = f.geometry();

  // tile_id in meters
  var p = geom.transform(REF_PROJ, 1);
  var xy = ee.List(p.coordinates());
  var x = ee.Number(xy.get(0));
  var y = ee.Number(xy.get(1));
  var ix = x.divide(GRID_SIZE).floor();
  var iy = y.divide(GRID_SIZE).floor();
  var tile_id = ix.format().cat('_').cat(iy.format());

  // lon/lat
  var ll = geom.transform('EPSG:4326', 1);
  var llxy = ee.List(ll.coordinates());
  var lon = ee.Number(llxy.get(0));
  var lat = ee.Number(llxy.get(1));

  return f.set({ tile_id: tile_id, Longitude: lon, Latitude: lat });
});

// QC filters
samples = samples
  .filter(ee.Filter.gte('valid_frac', VALID_FRAC_MIN))
  .filter(ee.Filter.gte('CHM_p95', CHM_P95_MIN))
  .filter(ee.Filter.lte('CHM_p95', CHM_P95_MAX));

// Optional thinning
if (THIN_FRAC < 1.0) {
  samples = samples.randomColumn('rand', 42).filter(ee.Filter.lt('rand', THIN_FRAC));
}

print('Pseudo samples (after QC):', samples.size());
print('Example rows:', samples.limit(5));

/**** =========================
 *  8) EXPORT CSV
 *  ========================= ****/

Export.table.toDrive({
  collection: samples,
  description: EXPORT_NAME,
  folder: DRIVE_FOLDER,
  fileFormat: 'CSV'
});
