// Script for collecting all farms GEE variables for a specified set/sequence of year(s)
// Paul Hegedus - 20191109
// This code runs in Google Earth Engine Code Editor.

//----------------------------------------------------------------
// Assets (COPY AND PASTE, UNCOMMENT, & IMPORT)
//----------------------------------------------------------------
//var gridmet = ee.ImageCollection("IDAHO_EPSCOR/GRIDMET"),
//ned = ee.Image("USGS/NED"),
//cdem = ee.ImageCollection("NRCan/CDEM"),
//srtm = ee.Image('USGS/SRTMGL1_003'),
//L7sr = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR"),
//s2toa = ee.ImageCollection("COPERNICUS/S2"),
//L5sr = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR"),
//L8sr = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR"),
//smap = ee.ImageCollection('NASA_USDA/HSL/SMAP_soil_moisture'),
//daymet = ee.ImageCollection('NASA/ORNL/DAYMET_V3'); // //s2sr = ee.ImageCollection("COPERNICUS/S2_SR"),
//olm_grtgroup = ee.Image("OpenLandMap/SOL/SOL_GRTGROUP_USDA-SOILTAX_C/v01"),
//olm_texture = ee.Image("OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02"),
//olm_bulkdensity = ee.Image("OpenLandMap/SOL/SOL_BULKDENS-FINEEARTH_USDA-4A1H_M/v02"),
//olm_claycontent = ee.Image("OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02"),
//olm_sandcontent = ee.Image("OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02"),
//olm_phw = ee.Image("OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02"),
//olm_watercontent = ee.Image("OpenLandMap/SOL/SOL_WATERCONTENT-33KPA_USDA-4B1C_M/v01"),
//olm_carboncontent = ee.Image("OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02");


// THIS VERSION COLLECTS DAYMET FOR USA FARMS b/c 1km RESOLUTION

// This code is for looping through a sequence of specified years for all farms 
// and downloading all GEE variables for a specified year. The code takes the user selected years
// and downloads all data for that year up until 03/30. It also downloads the full year's worth of 
// data for the previous year. This means that the least recent year of data in the database will be 
// 1 minus the least recent year specified in the user selection. The user can select a sequence of years
// in order, or a string of years to download data for. This code was developed for efficiently downloading
// previous years data from GEE. However, because the user can specify specific years, this code can be used as
// the annual run once per year (not ofpe_test/allFarms_singleYear data) and will download all the variables 
// needed for all the farms for the current year (up to mar) and previous year (full year) when 1 year is specified.
// In the database these are distinguished by a column labeled 'mar' or 'full' for up to march
// and full year, respectively. Note that water years are from november 1 to october 31. 
// User also has control over what variables to collect so not all vars are collected if not desired.

// Precipitation and growing degree day information is collected from GRIDMET by the Idaho EPSCoR project
// and is 4 km resolution with data available from 1979. The surfaces data is collected from the USGS
// National Elevation Dataset at a 10 m scale from prior to 1999. The vegetation index data is collected
// from different sources based on the date requested. All data 2018 or later are 10 m resolution filtered 
// top of atmosphere measures from the EU's Copernicus Sentinel 2. Prior to 2018 cy and 2016 py, data
// is collected from Landsat 5 if before 2012, Landsat 7 for 2012/2013, and Landsat 8 for after. 
// Landsat 8 is collected until present to gap fill. S2 provides NDRE and CLRE (should be CIRE) bands for calculating. 
// For Canadian farms the data is collected from NASA ORNL DaymetV3 at a 1km scale. These are estimates.
// Elevation data is from the National ELevaiton dataset at 10m scale for the USA and 
// from the Canadian Digital Elevation Model at a 20m scale for elevation and SRTM from NASA
// at a 30m scale for slope, aspect, tpi.

// NAMING CONVENTION
// The naming convention for these files are important for the import of the data to the database.
// The import code is relatively rigid because this GEE data is under our control and we can name
// and format how we like. The import code gathers information about the data from the name of the 
// file so the naming convention becomes important. It would be easiest if these were not changed.
// GENERAL FORMAT: farmname(sep="_")_measure_datasource_scale_year_(mar) <- only include if for current year, otherwise blank
// PREC & GDD FORMAT: farmname(sep="_")_scale_datasource_measure_(currYear/prevYear)_year <- currYear&prevYear are to indicate to mar or full
// e.g. “wood_carter_clre_20m_2013_mar” = 20 m resolution clre data from 2013 up until 03/30 for wood's fields near carter
// or   “merja_vaughn_4km_gdd_prevYr_2018” = 4km gdd data from 2018 for merja's fields near vaughn
// IMPORTANT: If you change the name of any variable in the filenames or if you use a new scale (i.e. not 4km,
// 10m, or 20m) you MUST EDIT CODE to change string patterns in functions. 


//----------------------------------------------------------------
// USER INPUTS
//----------------------------------------------------------------
// Select vars to collect 1=yes, 0=no
var precCY = 1; 
var precPY = 1; 
var gddCY = 1; 
var gddPY = 1; 
var aspect = 1; 
var slopeI = 1; 
var elev = 1; 
var tpiI = 1; 
var ndviCY = 1; 
var ndviPY = 1; 
var ndreCY = 1; 
var ndrePY = 1; 
var clreCY = 1; 
var clrePY = 1; 
var smapCY = 1;
var smapPY = 1;
var olmGrtgrp = 1;
var olmTexture = 1;
var olmBulkDens = 1;
var olmClayContent = 1;
var olmSandContent = 1;
var olmpHw = 1;
var olmWaterContent = 1;
var olmCarbonContent = 1;

// String or sequence of years to get data for
var years = ['2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'];
//var years = ['2015'];

// String of farms (assets that are ingested alread) to loop over (NEED to add others when you have them)
var farms = ['bailey_fortbenton','broyles_columbus','broyles_rapelje','loewen_steinbach','merja_cascade','merja_sunriver','merja_vaughn','norgaard_shonkin','oconnor_ekalaka','quinn_bigsandy','wood_carter','wood_loma','wood_virgelle','vandyke_gallatin']; //
//var farms = ['broyles_rapelje']; //
//****************************************************************
// START CODE
//****************************************************************
// start loop through years
//----------------------------------------------------------------
var j=0;
for(j; j < years.length; j++){
  var CY = years[j];
  var PY = CY - 1;
  // start for loop through farms
  //----------------------------------------------------------------
  var i=0; 
  for(i; i < farms.length; i++){
    //****************************************************************
    // set up assets
    //****************************************************************
    var table = ee.FeatureCollection("users/paulhegedus/" + farms[i] + "_bbox");
    var farm_name = table;
    var current_farm = farms[i];
    //Map.addLayer(farm_name,{}, current_farm + '_bbox', true);
    
    //********************************************************************
    // START CODE FOR FARMi YEARj
    //********************************************************************
    ///////////////////////////////////////////////////////////
    // CLIMATE VARIABLES (PREC, GDD)
    ///////////////////////////////////////////////////////////
    // Climate variables in North America colleceted from NASA ORNL
    // DaymetV3 dataset (1km)
    //----------------------------------------------------------------
    // Precipitation
    //----------------------------------------------------------------
    // PREC CY
    if(precCY==1){
      //use month-day of 11-01 and 03-30 for prec thru march of current grow yr, 
      var currYear = CY;
      var startDate = PY + '-11-01';
      var endDate =  CY + '-03-30';
      // select for precip data
      var precip = daymet.select('prcp')
        .filterDate(startDate, endDate)
        .sum();
      // Clip Precip images
      var mask = ee.FeatureCollection(farm_name);
      var mask_geometry = mask.geometry().bounds();
      var precip_sum_cy_clipped = precip.clip(mask);
      //export to google drive        
      Export.image.toDrive({
        image: precip_sum_cy_clipped,
        description: current_farm + '_1km_daymet_prec_currYr_' + currYear,
        folder:'GEE_Surfaces',
        crs: 'EPSG:4326', // saves it as long lat wgs84 (reprojected for spatial analysis)
        scale: 10,
        region: mask_geometry
      });
    } // END PREC CY
    //----------------------------------------------------------------
    // PREC PY
    if(precPY==1){
      //use 11-01 and 10-31 for prev grow yr
      var prevYear = PY;
      var startDate = PY-1 +'-11-01';
      var endDate =  PY + '-10-31';
      // select for precip data
      var precip = daymet.select('prcp')
        .filterDate(startDate, endDate)
        .sum();
      // Clip Precip images
      var mask = ee.FeatureCollection(farm_name);
      var mask_geometry = mask.geometry().bounds();
      var precip_sum_py_clipped = precip.clip(mask);
      //export to google drive
      Export.image.toDrive({
        image: precip_sum_py_clipped,
        description: current_farm + '_1km_daymet_prec_prevYr_' + prevYear,
        folder:'GEE_Surfaces',
        crs: 'EPSG:4326',
        scale: 10,
        region: mask_geometry
      });
    } // END PREC PY
    //----------------------------------------------------------------
    // Growing Degree Days
    //----------------------------------------------------------------
    // GDD CY
    if(gddCY==1){
      var grow_yr = CY;
      var startDate = CY + '-01-01';
      var endDate =  CY + '-03-30';
      // Calculate GDD Function
      var gddV3 = function(image){
        var days = image.expression('(tmax + tmin)/2 - 0', {
          'tmin': image.select('tmin'),
          'tmax': image.select('tmax')
        });
        return days.rename(['gdd']).updateMask(days.gt(0));
      };
      //use this block for gdd from jan 1 to mar 30 of current yr (current year)
      var gdd_cy = daymet 
        .filterDate(startDate, endDate)
        //.filterDate(grow_yr.toString(), (grow_yr + 1).toString())
        //.filter(ee.Filter.dayOfYear(1, 90))
        .map(gddV3)
        .sum();
      // Clip GDD images
      var mask = ee.FeatureCollection(farm_name);
      var mask_geometry = mask.geometry().bounds();
      var gdd_cy_clip =  gdd_cy.clip(mask); 
      //export to google drive
      Export.image.toDrive({
        image: gdd_cy_clip, //
        description: current_farm + '_1km_daymet_gdd_currYr_' + grow_yr, 
        folder:'GEE_Surfaces',
        crs: 'EPSG:4326',
        scale: 10,
        region: mask_geometry
      });
    } // END CY GDD
    //----------------------------------------------------------------
    // GDD PY
    if(gddPY==1){
      var grow_PrevYr = PY;
      var startDate = PY + '-01-01';
      var endDate =  PY + '-09-01';
      // Calculate GDD Function
      var gddV3 = function(image){
        var days = image.expression('(tmin + tmax)/2 - 0', {
          'tmin': image.select('tmin'),
          'tmax': image.select('tmax')
        });
        return days.rename(['gdd']).updateMask(days.gt(0));
      };
      // use this block for gdd from jan 1 thru jun 30 of prev grow yr (prev yr)
      var gdd_py = daymet
        .filterDate(startDate, endDate)
        //.filterDate(grow_PrevYr.toString(), (grow_PrevYr + 1).toString())
        //.filter(ee.Filter.dayOfYear(1, 180))
        .map(gddV3)
        .sum();
      // Clip GDD images
      var mask = ee.FeatureCollection(farm_name);
      var mask_geometry = mask.geometry().bounds();
      var gdd_py_clip = gdd_py.clip(mask); //
      //export to google drive
      Export.image.toDrive({
        image: gdd_py_clip, 
        description: current_farm + '_1km_daymet_gdd_prevYr_' + grow_PrevYr, 
        folder:'GEE_Surfaces',
        crs: 'EPSG:4326',
        scale: 10,
        region: mask_geometry
      });
    } // END PY GDD   
    //---------------------------------------------------------
    if(farms[i]!='loewen_steinbach'){ // IF NOT CANADA get GRIDMET
      //---------------------------------------------------------
      // Climate variables in the USA colleceted from UofIdaho
      // GRIDMET dataset (4km)
      //----------------------------------------------------------------
      // Precipitation
      //----------------------------------------------------------------
      // PREC CY
      if(precCY==1){
        //use month-day of 11-01 and 03-30 for prec thru march of current grow yr, 
        var currYear = CY;
        var startDate = PY + '-11-01';
        var endDate =  CY + '-03-30';
        // select for precip data
        var precip = gridmet.select('pr')
          .filterDate(startDate, endDate)
          .sum();
        // Clip Precip images
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        var precip_sum_cy_clipped = precip.clip(mask);
        //export to google drive        
        Export.image.toDrive({
          image: precip_sum_cy_clipped,
          description: current_farm + '_4km_gridmet_prec_currYr_' + currYear,
          folder:'GEE_Surfaces',
          crs: 'EPSG:4326', // saves it as long lat wgs84 (reprojected for spatial analysis)
          scale: 10,
          region: mask_geometry
        });
      } // END PREC CY
      //----------------------------------------------------------------
      // PREC PY
      if(precPY==1){
        //use 11-01 and 10-31 for prev grow yr
        var prevYear = PY;
        var startDate = PY-1 +'-11-01';
        var endDate =  PY + '-10-31';
        // select for precip data
        var precip = gridmet.select('pr')
          .filterDate(startDate, endDate)
          .sum();
        // Clip Precip images
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        var precip_sum_py_clipped = precip.clip(mask);
        //export to google drive
        Export.image.toDrive({
          image: precip_sum_py_clipped,
          description: current_farm + '_4km_gridmet_prec_prevYr_' + prevYear,
          folder:'GEE_Surfaces',
          crs: 'EPSG:4326',
          scale: 10,
          region: mask_geometry
        });
      } // END PREC PY
      //----------------------------------------------------------------
      // Growing Degree Days
      //----------------------------------------------------------------
      // GDD CY
      if(gddCY==1){
        var grow_yr = CY;
        var startDate = CY + '-01-01';
        var endDate =  CY + '-03-30';
        // Calculate GDD Function
        var gdd = function(image){
          var days = image.expression('(tmmn + tmmx)/2 - (0 + 273.15)', {
            'tmmn': image.select('tmmn'),
            'tmmx': image.select('tmmx')
          });
          return days.rename(['gdd']).updateMask(days.gt(0));
        };
        //use this block for gdd from jan 1 to mar 30 of current yr (current year)
        var gdd_cy = gridmet 
          .filterDate(startDate, endDate)
          //.filterDate(grow_yr.toString(), (grow_yr + 1).toString())
          //.filter(ee.Filter.dayOfYear(1, 90))
          .map(gdd)
          .sum();
        // Clip GDD images
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        var gdd_cy_clip =  gdd_cy.clip(mask); 
        //export to google drive
        Export.image.toDrive({
          image: gdd_cy_clip, //
          description: current_farm + '_4km_gridmet_gdd_currYr_' + grow_yr, 
          folder:'GEE_Surfaces',
          crs: 'EPSG:4326',
          scale: 10,
          region: mask_geometry
        });
      } // END CY GDD
      //----------------------------------------------------------------
      // GDD PY
      if(gddPY==1){
        var grow_PrevYr = PY;
        var startDate = PY + '-01-01';
        var endDate =  PY + '-09-01';
        // Calculate GDD Function
        var gdd = function(image){
          var days = image.expression('(tmmn + tmmx)/2 - (0 + 273.15)', {
            'tmmn': image.select('tmmn'),
            'tmmx': image.select('tmmx')
          });
          return days.rename(['gdd']).updateMask(days.gt(0));
        };
        // use this block for gdd from jan 1 thru jun 30 of prev grow yr (prev yr)
        var gdd_py = gridmet
          .filterDate(startDate, endDate)
          //.filterDate(grow_PrevYr.toString(), (grow_PrevYr + 1).toString())
          //.filter(ee.Filter.dayOfYear(1, 180))
          .map(gdd)
          .sum();
        // Clip GDD images
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        var gdd_py_clip = gdd_py.clip(mask); //
        //export to google drive
        Export.image.toDrive({
          image: gdd_py_clip, 
          description: current_farm + '_4km_gridmet_gdd_prevYr_' + grow_PrevYr, 
          folder:'GEE_Surfaces',
          crs: 'EPSG:4326',
          scale: 10,
          region: mask_geometry
        });
      } // END PY GDD
      //---------------------------------------------------------
    } // END IF USA
    //---------------------------------------------------------
    ///////////////////////////////////////////////////////////
    // VARIABLES FROM ELEVATION DATASETS (ELEV,SLOPE,ASPECT,TPI)
    ///////////////////////////////////////////////////////////
    if(farms[i]=='loewen_steinbach'){
      //---------------------------------------------------------
      // CANADIAN ELEVATION MODEL
      // 3/4 arc second resolution (~23 m)
      //---------------------------------------------------------
      var year = CY;
      // ELEVATION
      if(elev==1){
        var ced_ic = cdem.select('elevation');
        var ced = ced_ic.mean();
        var mask = farm_name;
        var mask_geometry = mask.geometry().bounds();
        var ced_clip = ced.clip(mask);
        //export to google drive
        Export.image.toDrive({
          image: ced_clip,  //
          description: current_farm + '_elev_cdem_20m_' + year, // 
          folder:'GEE_Surfaces',
          scale: 20,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      }
      // ASPECT_RAD
      if(aspect==1){
        // Define a function to convert from degrees to radians for aspect
        var radians = function(img) {
          return img.toFloat().multiply(Math.PI).divide(180);
        };
        var dem = srtm.select('elevation');
        var aspectDeg = ee.Terrain.aspect(dem);
        var aspect_rad = radians(aspectDeg.select('aspect')).float();
        var mask = farm_name;
        var mask_geometry = mask.geometry().bounds();
        var aspect_rad_clip = aspect_rad.clip(mask);
        //export to google drive
        Export.image.toDrive({
          image: aspect_rad_clip,  
          description: current_farm + '_aspect_rad_srtm_30m_' + year, 
          folder:'GEE_Surfaces',
          scale: 30,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      }
      // SLOPE
      if(slopeI==1){
        var dem = srtm.select('elevation');
        var slope = ee.Terrain.slope(dem).float();
        var mask = farm_name;
        var mask_geometry = mask.geometry().bounds();
        var slope_clip = slope.clip(mask);
        //export to google drive
        Export.image.toDrive({
          image: slope_clip,  //
          description: current_farm + '_slope_srtm_30m_' + year, // 
          folder:'GEE_Surfaces',
          scale: 30,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      }
      // TPI 30m
      if(tpiI==1){
        // calculate TPI
        var scale = 30; // scale in m
        var tpi_radius = 9*scale; // 3 * 30m = 270m radius.
        var tpi = srtm.subtract(srtm.focal_mean(tpi_radius, 'circle', 'meters')).add(0.5);
        var mask = farm_name;
        var mask_geometry = mask.geometry().bounds();
        var tpi_clip = tpi.clip(mask);
        //export to google drive
        Export.image.toDrive({
          image: tpi_clip,  // 
          description: current_farm + '_tpi_srtm_30m_' + year, //
          folder:'GEE_Surfaces',
          scale: 30,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      }
      //---------------------------------------------------------
    }else{ // END IF CANADIAN FARM
      //---------------------------------------------------------
      // USGS NATIONAL ELEVATION DATASET
      // downloads elevation data, calculates slope, aspect, tpi from USGS National 
      // Elevation Dataset (NED); 
      // @ 1/3 arc-second resolution (10m) clipped to farm boundary box
      //---------------------------------------------------------
      var year = CY;
      // ELEVATION
      if(elev==1){
        var mask = farm_name;
        var mask_geometry = mask.geometry().bounds();
        var dem_clip = ned.clip(mask);
        //export to google drive
        Export.image.toDrive({
          image: dem_clip,  //
          description: current_farm + '_elev_ned_10m_' + year, // 
          folder:'GEE_Surfaces',
          scale: 10,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      }
      // ASPECT_RAD
      if(aspect==1){
        // Define a function to convert from degrees to radians for aspect
        var radians = function(img) {
          return img.toFloat().multiply(Math.PI).divide(180);
        };
        var terrain = ee.Algorithms.Terrain(ned); // Calculates slope, aspect, and a simple hillshade from a terrain DEM.
        var aspect_rad = radians(terrain.select('aspect')).float();
        var mask = farm_name;
        var mask_geometry = mask.geometry().bounds();
        var aspect_rad_clip = aspect_rad.clip(mask);
        //export to google drive
        Export.image.toDrive({
          image: aspect_rad_clip,  
          description: current_farm + '_aspect_rad_ned_10m_' + year, 
          folder:'GEE_Surfaces',
          scale: 10,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      }
      // SLOPE
      if(slopeI==1){
        var terrain = ee.Algorithms.Terrain(ned); // Calculates slope, aspect, and a simple hillshade from a terrain DEM.
        var slope = terrain.select('slope').float();
        var mask = farm_name;
        var mask_geometry = mask.geometry().bounds();
        var slope_clip = slope.clip(mask);
        //export to google drive
        Export.image.toDrive({
          image: slope_clip,  //
          description: current_farm + '_slope_ned_10m_' + year, // 
          folder:'GEE_Surfaces',
          scale: 10,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      }
      // TPI 30m
      if(tpiI==1){
        // calculate TPI
        var scale = 10; // scale in m
        var tpi_radius = 3*scale; // 3 * 10m = 30m res.
        var tpi = ned.subtract(ned.focal_mean(tpi_radius, 'circle', 'meters')).add(0.5);
        var mask = farm_name;
        var mask_geometry = mask.geometry().bounds();
        var tpi_clip = tpi.clip(mask);
        //export to google drive
        Export.image.toDrive({
          image: tpi_clip,  // 
          description: current_farm + '_tpi_ned_10m_' + year, //
          folder:'GEE_Surfaces',
          scale: 10,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      }
      //---------------------------------------------------------
    } // END IF USA FARM
    //---------------------------------------------------------
    ///////////////////////////////////////////////////////////
    // SMAP - Soil Moisture Active Passive (20km )
    ///////////////////////////////////////////////////////////
    //----------------------------------------------------
    if(smapCY==1){
      var year = CY; //  
      if(CY >= 2016){
        var currYear = CY;
        var startDate = PY + '-11-01';
        var endDate =  CY + '-03-30';
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        // SURFACE SOIL MOISTURE (mm) 20km @45N
        var ssmCY = smap.select('ssm') //.select('ssm')
          .filterDate(startDate, endDate)
          .filterBounds(mask_geometry)
          .map(function(image){return image.clip(mask_geometry)}) ;
        var ssmCY = ssmCY.reduce(ee.Reducer.median());
        // SUBSURFACE SOIL MOISTURE (mm) 20km @45N
        var susmCY = smap.select('susm')
          .filterDate(startDate, endDate)
          .filterBounds(mask_geometry)
          .map(function(image){return image.clip(mask_geometry)}) ;
        var susmCY = susmCY.reduce(ee.Reducer.median());
        // export to google drive
        Export.image.toDrive({
          image: ssmCY,
          description: current_farm + '_ssm_smap_20km'  + '_' + year + '_mar',
          folder:'GEE_Surfaces',
          scale: 10,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: susmCY,
          description: current_farm + '_susm_smap_20km'  + '_' + year + '_mar',
          folder:'GEE_Surfaces',
          scale: 10,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // END IF YEAR >= 2016
      //---------------------------------------------------------
    } // END CY SMAP
    //----------------------------------------------------
    if(smapPY==1){
      var year = PY; //  
      if(CY >= 2016){
        var prevYear = PY;
        var startDate = PY-1 +'-11-01';
        var endDate =  PY + '-10-31';
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        // SURFACE SOIL MOISTURE (mm) 20km @45N
        var ssmPY = smap.select('ssm')
          .filterDate(startDate, endDate)
          .filterBounds(mask_geometry)
          .map(function(image){return image.clip(mask_geometry)}) ;
        var ssmPY = ssmPY.reduce(ee.Reducer.median());
        // SUBSURFACE SOIL MOISTURE (mm) 20km @45N
        var susmPY = smap.select('susm')
          .filterDate(startDate, endDate)
          .filterBounds(mask_geometry)
          .map(function(image){return image.clip(mask_geometry)}) ;
        var susmPY = susmPY.reduce(ee.Reducer.median());
        // export to google drive
        Export.image.toDrive({
          image: ssmPY,
          description: current_farm + '_ssm_smap_20km'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 10,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: susmPY,
          description: current_farm + '_susm_smap_20km'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 10,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // END IF YEAR >= 2016
      //---------------------------------------------------------
    } // END PY SMAP
    //----------------------------------------------------
    ///////////////////////////////////////////////////////////
    // NDVI 
    ///////////////////////////////////////////////////////////
    //---------------------------------------------------------------
    if(ndviCY==1){
      var year = CY; 
      var start = ee.Date(year+'-01-01');
      var end = ee.Date(year+'-03-30'); //
      // based on year decide what data source to use
      if(CY >= 2016){
        // calculates ndvi @ 10m resolution from Sentinel-2 
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // Sentinel 2 (10 m)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // s2_filter: filters all avail sentinel 2 (s2) images: limit spatial extent 
        // to the boundary and filter date range by start and end vars defined above
        // Function to mask clouds using the Sentinel-2 QA band.
        var maskS2clouds = function(image) {
          var qa = image.select('QA60');
          // Bits 10 and 11 are clouds and cirrus, respectively.
          var cloudBitMask = 1 << 10;
          var cirrusBitMask = 1 << 11;
          // Both flags should be set to zero, indicating clear conditions.
          var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
                     qa.bitwiseAnd(cirrusBitMask).eq(0));
          // Return the masked and scaled data, without the QA bands.
          return image.updateMask(mask).divide(10000)
              .select("B.*")
              .copyProperties(image, ["system:time_start"]);
        };
        var s2toa_filterCY = s2toa
          .filterBounds(farm_name)
          .filterDate(start, end)
          // Pre-filter to get less cloudy granules.
          .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
          .map(maskS2clouds);
        // compute ndvi using red band 4 and NIR band 8a; calculated as
        // : (b8-b4)/(b8+b4)
        var ndviTOACY = s2toa_filterCY.map(function(img){
          var red = ee.Image(img.select('B4'));
          var nir = ee.Image(img.select('B8'));
          return (nir.subtract(red)).divide(nir.add(red)).rename('ndvi');
        });
        // Make a "greenest" pixel composite.
        var ndvi_95TOACY = ndviTOACY.qualityMosaic('ndvi');
        // clip images
        var mask = ee.FeatureCollection(farm_name);
        var ndvi_95_clipTOACY = ndvi_95TOACY.clip(mask);
        var mask_geometry = mask.geometry().bounds();
        // export to google drive
        Export.image.toDrive({
          image: ndvi_95_clipTOACY,
          description: current_farm + '_ndvi_S2_10m'  + '_' + year + '_mar',
          folder:'GEE_Surfaces',
          scale: 10,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // Landsat 8 (10 m)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        var maskL8sr = function(image) {
            // Bits 3 and 5 are cloud shadow and cloud, respectively.
            var cloudShadowBitMask = (1 << 3);
            var cloudsBitMask = (1 << 5);
            // Get the pixel QA band.
            var qa = image.select('pixel_qa');
            // Both flags should be set to zero, indicating clear conditions.
            var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                      .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
            return image.updateMask(mask);
          };
          var L8SRCY_filter = L8sr
            .filterBounds(farm_name)
            .filterDate(start, end)
            .map(maskL8sr);
          var ndviSRCYL8 = L8SRCY_filter.map(function(img){
            var red = ee.Image(img.select('B4'));
            var nir = ee.Image(img.select('B5'));
            return (nir.subtract(red)).divide(nir.add(red)).rename('ndvi');
          });
          // Make a "greenest" pixel composite.
          var ndvi_95SRCYL8 = ndviSRCYL8.qualityMosaic('ndvi');
          // clip images
          var mask = ee.FeatureCollection(farm_name);
          var ndvi_95_clipSRCYL8 = ndvi_95SRCYL8.clip(mask);
          var mask_geometry = mask.geometry().bounds();
          // export to google drive
          Export.image.toDrive({
            image: ndvi_95_clipSRCYL8,
            description: current_farm + '_ndvi_L8_30m'  + '_' + year + '_mar',
            folder:'GEE_Surfaces',
            scale: 30,
            crs: 'EPSG:4326',
            region: mask_geometry
          });
      }else{ // END IF CY >= 2016
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // LANDSAT 5 (1999py to 2012cy)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if(CY <= 2012){ 
          // calculates ndvi @ 30m resolution from Landsat 5 surface reflectance
          var cloudMaskL457 = function(image) {
            var qa = image.select('pixel_qa');
            // If the cloud bit (5) is set and the cloud confidence (7) is high
            // or the cloud shadow bit is set (3), then it's a bad pixel.
            var cloud = qa.bitwiseAnd(1 << 5)
                      .and(qa.bitwiseAnd(1 << 7))
                      .or(qa.bitwiseAnd(1 << 3));
            // Remove edge pixels that don't occur in all bands
            var mask2 = image.mask().reduce(ee.Reducer.min());
            return image.updateMask(cloud.not()).updateMask(mask2);
          };
          var L5SRCY_filter = L5sr
            .filterBounds(farm_name)
            .filterDate(start, end)
            .map(cloudMaskL457);
          var ndviSRCYL5 = L5SRCY_filter.map(function(img){
            var red = ee.Image(img.select('B3'));
            var nir = ee.Image(img.select('B4'));
            return (nir.subtract(red)).divide(nir.add(red)).rename('ndvi');
          });
          // Make a "greenest" pixel composite.
          var ndvi_95SRCYL5 = ndviSRCYL5.qualityMosaic('ndvi');
          // clip images
          var mask = ee.FeatureCollection(farm_name);
          var ndvi_95_clipSRCYL5 = ndvi_95SRCYL5.clip(mask);
          var mask_geometry = mask.geometry().bounds();
          // export to google drive
          Export.image.toDrive({
            image: ndvi_95_clipSRCYL5,
            description: current_farm + '_ndvi_L5_30m'  + '_' + year + '_mar',
            folder:'GEE_Surfaces',
            scale: 30,
            crs: 'EPSG:4326',
            region: mask_geometry
          });
        } // END L5 IF CY <=2012
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // LANDSAT 7 (2012py & 2013cy)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if(CY == 2013){ 
          // calculates ndvi @ 30m resolution from Landsat 7 surface reflectance
          var cloudMaskL457 = function(image) {
            var qa = image.select('pixel_qa');
            // If the cloud bit (5) is set and the cloud confidence (7) is high
            // or the cloud shadow bit is set (3), then it's a bad pixel.
            var cloud = qa.bitwiseAnd(1 << 5)
                      .and(qa.bitwiseAnd(1 << 7))
                      .or(qa.bitwiseAnd(1 << 3));
            // Remove edge pixels that don't occur in all bands
            var mask2 = image.mask().reduce(ee.Reducer.min());
            return image.updateMask(cloud.not()).updateMask(mask2);
          };
          var L7SRCY_filter = L7sr
            .filterBounds(farm_name)
            .filterDate(start, end)
            .map(cloudMaskL457);
          var ndviSRCYL7 = L7SRCY_filter.map(function(img){
            var red = ee.Image(img.select('B3'));
            var nir = ee.Image(img.select('B4'));
            return (nir.subtract(red)).divide(nir.add(red)).rename('ndvi');
          });
          // Make a "greenest" pixel composite.
          var ndvi_95SRCYL7 = ndviSRCYL7.qualityMosaic('ndvi');
          // clip images
          var mask = ee.FeatureCollection(farm_name);
          var ndvi_95_clipSRCYL7 = ndvi_95SRCYL7.clip(mask);
          var mask_geometry = mask.geometry().bounds();
          // export to google drive
          Export.image.toDrive({
            image: ndvi_95_clipSRCYL7,
            description: current_farm + '_ndvi_L7_30m'  + '_' + year + '_mar',
            folder:'GEE_Surfaces',
            scale: 30,
            crs: 'EPSG:4326',
            region: mask_geometry
          });
        } // END L7 IF CY == 2013
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // LANDSAT 8 (2013py (april - dec) to 2018cy or 2016py)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if(CY >= 2014){ 
          var maskL8sr = function(image) {
            // Bits 3 and 5 are cloud shadow and cloud, respectively.
            var cloudShadowBitMask = (1 << 3);
            var cloudsBitMask = (1 << 5);
            // Get the pixel QA band.
            var qa = image.select('pixel_qa');
            // Both flags should be set to zero, indicating clear conditions.
            var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                      .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
            return image.updateMask(mask);
          };
          var L8SRCY_filter = L8sr
            .filterBounds(farm_name)
            .filterDate(start, end)
            .map(maskL8sr);
          var ndviSRCYL8 = L8SRCY_filter.map(function(img){
            var red = ee.Image(img.select('B4'));
            var nir = ee.Image(img.select('B5'));
            return (nir.subtract(red)).divide(nir.add(red)).rename('ndvi');
          });
          // Make a "greenest" pixel composite.
          var ndvi_95SRCYL8 = ndviSRCYL8.qualityMosaic('ndvi');
          // clip images
          var mask = ee.FeatureCollection(farm_name);
          var ndvi_95_clipSRCYL8 = ndvi_95SRCYL8.clip(mask);
          var mask_geometry = mask.geometry().bounds();
          // export to google drive
          Export.image.toDrive({
            image: ndvi_95_clipSRCYL8,
            description: current_farm + '_ndvi_L8_30m'  + '_' + year + '_mar',
            folder:'GEE_Surfaces',
            scale: 30,
            crs: 'EPSG:4326',
            region: mask_geometry
          });
        } // END L8 IF CY >= 2014 AND  CY <= 2017
      } // END IF NOT >= 2018 S2
    } // end if current year ndvi
    //---------------------------------------------------------------
    // Previous year (full year)
    //---------------------------------------------------------------
    if(ndviPY==1){
      var year = PY; //  
      var start = ee.Date(year+'-01-01');
      var end = ee.Date(year+'-12-31'); //
      if(PY >= 2016){
        // calculates ndvi @ 10m resolution from Sentinel-2 
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // Sentinel 2 (10m)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // s2_filter: filters all avail sentinel 2 (s2) images: limit spatial extent 
        // to the boundary and filter date range by start and end vars defined above
        // Function to mask clouds using the Sentinel-2 QA band.
        var maskS2clouds = function(image) {
          var qa = image.select('QA60');
          // Bits 10 and 11 are clouds and cirrus, respectively.
          var cloudBitMask = 1 << 10;
          var cirrusBitMask = 1 << 11;
          // Both flags should be set to zero, indicating clear conditions.
          var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
                     qa.bitwiseAnd(cirrusBitMask).eq(0));
          // Return the masked and scaled data, without the QA bands.
          return image.updateMask(mask).divide(10000)
              .select("B.*")
              .copyProperties(image, ["system:time_start"]);
        };
        var s2toa_filterPY = s2toa
          .filterBounds(farm_name)
          .filterDate(start, end)
          // Pre-filter to get less cloudy granules.
          .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
          .map(maskS2clouds);
        // compute ndvi using red band 4 and NIR band 8a; calculated as
        // : (b8-b4)/(b8+b4)
        var ndviTOAPY = s2toa_filterPY.map(function(img){
          var red = ee.Image(img.select('B4'));
          var nir = ee.Image(img.select('B8'));
          return (nir.subtract(red)).divide(nir.add(red)).rename('ndvi');
        });
        // Make a "greenest" pixel composite.
        var ndvi_95TOAPY = ndviTOAPY.qualityMosaic('ndvi');
        // clip images
        var mask = ee.FeatureCollection(farm_name);
        var ndvi_95_clipTOAPY = ndvi_95TOAPY.clip(mask);
        var mask_geometry = mask.geometry().bounds();
        //export to google drive
        Export.image.toDrive({
          image: ndvi_95_clipTOAPY,
          description: current_farm + '_ndvi_S2_10m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 10,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // Landsat 8 (30m)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        var maskL8sr = function(image) {
          // Bits 3 and 5 are cloud shadow and cloud, respectively.
          var cloudShadowBitMask = (1 << 3);
          var cloudsBitMask = (1 << 5);
          // Get the pixel QA band.
          var qa = image.select('pixel_qa');
          // Both flags should be set to zero, indicating clear conditions.
          var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                  .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
          return image.updateMask(mask);
        };
        var L8SRPY_filter = L8sr
          .filterBounds(farm_name)
          .filterDate(start, end)
          .map(maskL8sr);
        var ndviSRPYfilterL8 = L8SRPY_filter.map(function(img){
          var red = ee.Image(img.select('B4'));
          var nir = ee.Image(img.select('B5'));
          return (nir.subtract(red)).divide(nir.add(red)).rename('ndvi');
        });
        // Make a "greenest" pixel composite.
        var ndvi_95SRPYfilterL8 = ndviSRPYfilterL8.qualityMosaic('ndvi');
        // clip images
        var mask = ee.FeatureCollection(farm_name);
        var ndvi_95_clipSRPYfilterL8 = ndvi_95SRPYfilterL8.clip(mask);
        var mask_geometry = mask.geometry().bounds();
        // export to google drive
        Export.image.toDrive({
          image: ndvi_95_clipSRPYfilterL8,
          description: current_farm + '_ndvi_L8_30m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 30,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      }else{ // END S2 PY >= 2016
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // LANDSAT 5 (1999py & 2012cy)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if(PY <= 2011){ 
          // calculates ndvi @ 30m resolution from Landsat 5 surface reflectance
          var cloudMaskL457 = function(image) {
            var qa = image.select('pixel_qa');
            // If the cloud bit (5) is set and the cloud confidence (7) is high
            // or the cloud shadow bit is set (3), then it's a bad pixel.
            var cloud = qa.bitwiseAnd(1 << 5)
                      .and(qa.bitwiseAnd(1 << 7))
                      .or(qa.bitwiseAnd(1 << 3));
            // Remove edge pixels that don't occur in all bands
            var mask2 = image.mask().reduce(ee.Reducer.min());
            return image.updateMask(cloud.not()).updateMask(mask2);
          };
          var L5SRPY_filter = L5sr
            .filterBounds(farm_name)
            .filterDate(start, end)
            .map(cloudMaskL457);
          var ndviSRPYfilterL5 = L5SRPY_filter.map(function(img){
            var red = ee.Image(img.select('B3'));
            var nir = ee.Image(img.select('B4'));
            return (nir.subtract(red)).divide(nir.add(red)).rename('ndvi');
          });
          // Make a "greenest" pixel composite.
          var ndvi_95SRPYfilterL5 = ndviSRPYfilterL5.qualityMosaic('ndvi');
          // clip images
          var mask = ee.FeatureCollection(farm_name);
          var ndvi_95_clipSRPYfilterL5 = ndvi_95SRPYfilterL5.clip(mask);
          var mask_geometry = mask.geometry().bounds();
          // export to google drive
          Export.image.toDrive({
            image: ndvi_95_clipSRPYfilterL5,
            description: current_farm + '_ndvi_L5_30m'  + '_' + year,
            folder:'GEE_Surfaces',
            scale: 30,
            crs: 'EPSG:4326',
            region: mask_geometry
          });
        } // END L5 FOR 1999 to 2011 PY
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // LANDSAT 7 (2012py & 2013cy)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if(PY == 2012){ 
          // calculates ndvi @ 30m resolution from Landsat 7 surface reflectance
          var cloudMaskL457 = function(image) {
            var qa = image.select('pixel_qa');
            // If the cloud bit (5) is set and the cloud confidence (7) is high
            // or the cloud shadow bit is set (3), then it's a bad pixel.
            var cloud = qa.bitwiseAnd(1 << 5)
                      .and(qa.bitwiseAnd(1 << 7))
                      .or(qa.bitwiseAnd(1 << 3));
            // Remove edge pixels that don't occur in all bands
            var mask2 = image.mask().reduce(ee.Reducer.min());
            return image.updateMask(cloud.not()).updateMask(mask2);
          };
          var L7SRPY_filter = L7sr
            .filterBounds(farm_name)
            .filterDate(start, end)
            .map(cloudMaskL457);
          var ndviSRPYfilterL7 = L7SRPY_filter.map(function(img){
            var red = ee.Image(img.select('B3'));
            var nir = ee.Image(img.select('B4'));
            return (nir.subtract(red)).divide(nir.add(red)).rename('ndvi');
          });
          // Make a "greenest" pixel composite.
          var ndvi_95SRPYfilterL7 = ndviSRPYfilterL7.qualityMosaic('ndvi');
          // clip images
          var mask = ee.FeatureCollection(farm_name);
          var ndvi_95_clipSRPYfilterL7 = ndvi_95SRPYfilterL7.clip(mask);
          var mask_geometry = mask.geometry().bounds();
          // export to google drive
          Export.image.toDrive({
            image: ndvi_95_clipSRPYfilterL7,
            description: current_farm + '_ndvi_L7_30m'  + '_' + year,
            folder:'GEE_Surfaces',
            scale: 30,
            crs: 'EPSG:4326',
            region: mask_geometry
          });
        } // END L7 FOR 2012 PY
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // LANDSAT 8 (2013py (april - dec) to 2018cy or 2016py)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if(PY >= 2013){ 
          if(PY==2013){
            //CHANGE DATES TO APRIL TO DEC
            var start13 = ee.Date(year+'-04-01');
            var end13 = ee.Date(year+'-12-31'); //
            var maskL8sr = function(image) {
              // Bits 3 and 5 are cloud shadow and cloud, respectively.
              var cloudShadowBitMask = (1 << 3);
              var cloudsBitMask = (1 << 5);
              // Get the pixel QA band.
              var qa = image.select('pixel_qa');
              // Both flags should be set to zero, indicating clear conditions.
              var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                      .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
              return image.updateMask(mask);
            };
            var L8SRPY_filter = L8sr
              .filterBounds(farm_name)
              .filterDate(start13, end13)
              .map(maskL8sr);
            var ndviSRPYfilterL8 = L8SRPY_filter.map(function(img){
              var red = ee.Image(img.select('B4'));
              var nir = ee.Image(img.select('B5'));
              return (nir.subtract(red)).divide(nir.add(red)).rename('ndvi');
            });
            // Make a "greenest" pixel composite.
            var ndvi_95SRPYfilterL8 = ndviSRPYfilterL8.qualityMosaic('ndvi');
            // clip images
            var mask = ee.FeatureCollection(farm_name);
            var ndvi_95_clipSRPYfilterL8 = ndvi_95SRPYfilterL8.clip(mask);
            var mask_geometry = mask.geometry().bounds();
            // export to google drive
            Export.image.toDrive({
              image: ndvi_95_clipSRPYfilterL8,
              description: current_farm + '_ndvi_L8_30m'  + '_' + year,
              folder:'GEE_Surfaces',
              scale: 30,
              crs: 'EPSG:4326',
              region: mask_geometry
            });
          }else{ // end if PY == 2013
            // USE NORMAL DATES (JAN TO DEC)
            var maskL8sr = function(image) {
              // Bits 3 and 5 are cloud shadow and cloud, respectively.
              var cloudShadowBitMask = (1 << 3);
              var cloudsBitMask = (1 << 5);
              // Get the pixel QA band.
              var qa = image.select('pixel_qa');
              // Both flags should be set to zero, indicating clear conditions.
              var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                      .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
              return image.updateMask(mask);
            };
            var L8SRPY_filter = L8sr
              .filterBounds(farm_name)
              .filterDate(start, end)
              .map(maskL8sr);
            var ndviSRPYfilterL8 = L8SRPY_filter.map(function(img){
              var red = ee.Image(img.select('B4'));
              var nir = ee.Image(img.select('B5'));
              return (nir.subtract(red)).divide(nir.add(red)).rename('ndvi');
            });
            // Make a "greenest" pixel composite.
            var ndvi_95SRPYfilterL8 = ndviSRPYfilterL8.qualityMosaic('ndvi');
            // clip images
            var mask = ee.FeatureCollection(farm_name);
            var ndvi_95_clipSRPYfilterL8 = ndvi_95SRPYfilterL8.clip(mask);
            var mask_geometry = mask.geometry().bounds();
            // export to google drive
            Export.image.toDrive({
              image: ndvi_95_clipSRPYfilterL8,
              description: current_farm + '_ndvi_L8_30m'  + '_' + year,
              folder:'GEE_Surfaces',
              scale: 30,
              crs: 'EPSG:4326',
              region: mask_geometry
            });
          } // end if not py == 2013 (apr to dec L8)
        } // END L8 FOR 2013 to 2015 PY
      } // END S2 ELSE PY >= 2016 
    } // END PY NDVI
    //---------------------------------------------------------------
    ///////////////////////////////////////////////////////////
    // NDRE 
    ///////////////////////////////////////////////////////////
    //---------------------------------------------------------------
    if(ndreCY==1){
      var year = CY; 
      var start = ee.Date(year+'-01-01');
      var end = ee.Date(year+'-03-30'); //
      // based on year decide what data source to use
      if(CY >= 2016){
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // Sentinel 2 (10m)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // calculates ndre @ 10m resolution from Sentinel-2 
        // s2_filter: filters all avail sentinel 2 (s2) images: limit spatial extent 
        // to the boundary and filter date range by start and end vars defined above
        // Function to mask clouds using the Sentinel-2 QA band.
        var maskS2clouds = function(image) {
          var qa = image.select('QA60');
          // Bits 10 and 11 are clouds and cirrus, respectively.
          var cloudBitMask = 1 << 10;
          var cirrusBitMask = 1 << 11;
          // Both flags should be set to zero, indicating clear conditions.
          var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
                     qa.bitwiseAnd(cirrusBitMask).eq(0));
          // Return the masked and scaled data, without the QA bands.
          return image.updateMask(mask).divide(10000)
              .select("B.*")
              .copyProperties(image, ["system:time_start"]);
        };
        var s2toa_filterCY = s2toa
          .filterBounds(farm_name)
          .filterDate(start, end)
          // Pre-filter to get less cloudy granules.
          .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
          .map(maskS2clouds);
        // compute ndre using red edge band 1 and red edge band 2; calculated as
        // : (b6-b5)/(b6+b5)
        var ndreTOACY = s2toa_filterCY.map(function(img){
          var re1 = ee.Image(img.select('B5'));
          var re2 = ee.Image(img.select('B6'));
          return (re1.subtract(re2)).divide(re1.add(re2)).rename('ndre');
        });
        // Make a "greenest" pixel composite.
        var ndre_95TOACY = ndreTOACY.qualityMosaic('ndre');
        // clip images
        var mask = ee.FeatureCollection(farm_name);
        var ndre_95_clipTOACY = ndre_95TOACY.clip(mask);
        var mask_geometry = mask.geometry().bounds();
        // export to google drive
        Export.image.toDrive({
          image: ndre_95_clipTOACY,
          description: current_farm + '_ndre_S2_20m'  + '_' + year + '_mar',
          folder:'GEE_Surfaces',
          scale: 20,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // end if year >= 2018
    } // end if current year ndre
    //---------------------------------------------------------------
    // Previous year (full year)
    //---------------------------------------------------------------
    if(ndrePY==1){
      var year = PY; //  
      var start = ee.Date(year+'-01-01');
      var end = ee.Date(year+'-12-31'); //
      if(PY >= 2016){
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // Sentinel 2 (10m)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // calculates ndre @ 10m resolution from Sentinel-2 
        // s2_filter: filters all avail sentinel 2 (s2) images: limit spatial extent 
        // to the boundary and filter date range by start and end vars defined above
        ///// PY TOA S2 w/filter /////
        // Function to mask clouds using the Sentinel-2 QA band.
        var maskS2clouds = function(image) {
          var qa = image.select('QA60');
          // Bits 10 and 11 are clouds and cirrus, respectively.
          var cloudBitMask = 1 << 10;
          var cirrusBitMask = 1 << 11;
          // Both flags should be set to zero, indicating clear conditions.
          var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
                     qa.bitwiseAnd(cirrusBitMask).eq(0));
          // Return the masked and scaled data, without the QA bands.
          return image.updateMask(mask).divide(10000)
              .select("B.*")
              .copyProperties(image, ["system:time_start"]);
        };
        var s2toa_filterPY = s2toa
          .filterBounds(farm_name)
          .filterDate(start, end)
          // Pre-filter to get less cloudy granules.
          .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
          .map(maskS2clouds);
        // compute ndre using red band 4 and NIR band 8a; calculated as
        // : (b8-b4)/(b8+b4)
        var ndreTOAPY = s2toa_filterPY.map(function(img){
          var re1 = ee.Image(img.select('B5'));
          var re2 = ee.Image(img.select('B6'));
          return (re1.subtract(re2)).divide(re1.add(re2)).rename('ndre');
        });
        // Make a "greenest" pixel composite.
        var ndre_95TOAPY = ndreTOAPY.qualityMosaic('ndre');
        // clip images
        var mask = ee.FeatureCollection(farm_name);
        var ndre_95_clipTOAPY = ndre_95TOAPY.clip(mask);
        var mask_geometry = mask.geometry().bounds();
        // export to google drive
        Export.image.toDrive({
          image: ndre_95_clipTOAPY,
          description: current_farm + '_ndre_S2_20m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 20,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // end if py >= 2016
    } // !!! END PY NDRE !!!  
    //---------------------------------------------------------------
    ///////////////////////////////////////////////////////////
    // CIRE 
    ///////////////////////////////////////////////////////////
    //---------------------------------------------------------------
    if(clreCY==1){
      var year = CY; 
      var start = ee.Date(year+'-01-01');
      var end = ee.Date(year+'-03-30'); //
      // based on year decide what data source to use
      if(CY >= 2016){
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // Sentinel 2 (10m)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // calculates clre @ 10m resolution from Sentinel-2 
        // s2_filter: filters all avail sentinel 2 (s2) images: limit spatial extent 
        // to the boundary and filter date range by start and end vars defined above
        // Function to mask clouds using the Sentinel-2 QA band.
        var maskS2clouds = function(image) {
          var qa = image.select('QA60');
          // Bits 10 and 11 are clouds and cirrus, respectively.
          var cloudBitMask = 1 << 10;
          var cirrusBitMask = 1 << 11;
          // Both flags should be set to zero, indicating clear conditions.
          var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
                     qa.bitwiseAnd(cirrusBitMask).eq(0));
          // Return the masked and scaled data, without the QA bands.
          return image.updateMask(mask).divide(10000)
              .select("B.*")
              .copyProperties(image, ["system:time_start"]);
        };
        var s2toa_filterCY = s2toa
          .filterBounds(farm_name)
          .filterDate(start, end)
          // Pre-filter to get less cloudy granules.
          .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
          .map(maskS2clouds);
        // calculate chlorophyll red edge index using bands 7 and 5; formula from 
        // clevers and gitleson 2012: (b7 / b5) - 1
        var clreTOACY = s2toa_filterCY.map(function(img){
          var b7 = ee.Image(img.select('B7'));
          var b5 = ee.Image(img.select('B5'));
          return (b7.divide(b5)).subtract(1).rename('clre');
        });
        // Make a "greenest" pixel composite.
        var clre_95TOACY = clreTOACY.qualityMosaic('clre');
        // clip images
        var mask = ee.FeatureCollection(farm_name);
        var clre_95_clipTOACY = clre_95TOACY.clip(mask);
        var mask_geometry = mask.geometry().bounds();
        // export to google drive
        Export.image.toDrive({
          image: clre_95_clipTOACY,
          description: current_farm + '_clre_S2_20m'  + '_' + year + '_mar',
          folder:'GEE_Surfaces',
          scale: 20,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // end if year >= 2016
    } // end if current year clre
    //---------------------------------------------------------------
    // Previous year (full year)
    //---------------------------------------------------------------
    if(clrePY==1){
      var year = PY; //  
      var start = ee.Date(year+'-01-01');
      var end = ee.Date(year+'-12-31'); //
      if(PY >= 2016){
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // Sentinel 2 (10m)
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // calculates clre @ 10m resolution from Sentinel-2 
        // s2_filter: filters all avail sentinel 2 (s2) images: limit spatial extent 
        // to the boundary and filter date range by start and end vars defined above
        // Function to mask clouds using the Sentinel-2 QA band.
        var maskS2clouds = function(image) {
          var qa = image.select('QA60');
          // Bits 10 and 11 are clouds and cirrus, respectively.
          var cloudBitMask = 1 << 10;
          var cirrusBitMask = 1 << 11;
          // Both flags should be set to zero, indicating clear conditions.
          var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
                     qa.bitwiseAnd(cirrusBitMask).eq(0));
          // Return the masked and scaled data, without the QA bands.
          return image.updateMask(mask).divide(10000)
              .select("B.*")
              .copyProperties(image, ["system:time_start"]);
        };
        var s2toa_filterPY = s2toa
          .filterBounds(farm_name)
          .filterDate(start, end)
          // Pre-filter to get less cloudy granules.
          .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
          .map(maskS2clouds);
        // calculate chlorophyll red edge index using bands 7 and 5; formula from 
        // clevers and gitleson 2012: (b7 / b5) - 1
        var clreTOAPY = s2toa_filterPY.map(function(img){
          var b7 = ee.Image(img.select('B7'));
          var b5 = ee.Image(img.select('B5'));
          return (b7.divide(b5)).subtract(1).rename('clre');
        });
        // Make a "greenest" pixel composite.
        var clre_95TOAPY = clreTOAPY.qualityMosaic('clre');
        // clip images
        var mask = ee.FeatureCollection(farm_name);
        var clre_95_clipTOAPY = clre_95TOAPY.clip(mask);
        var mask_geometry = mask.geometry().bounds();
        // export to google drive
        Export.image.toDrive({
          image: clre_95_clipTOAPY,
          description: current_farm + '_clre_S2_20m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 20,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // end if py >= 2016
    } // END PY CIRE
    //---------------------------------------------------------------
    ///////////////////////////////////////////////////////////
    // OPEN LAND MAP - GREAT GROUPS - (250m)
    ///////////////////////////////////////////////////////////
    //----------------------------------------------------
    if(olmGrtgrp==1){
      var year = CY; //  
      if(CY <= 2018){
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        Export.image.toDrive({
          image: olm_grtgroup,
          description: current_farm + '_olm_grtgroup_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // END IF YEAR <= 2018
    } // END GRTGROUP
    
    ///////////////////////////////////////////////////////////
    // OPEN LAND MAP - SOIL TEXTURE - (250m)
    ///////////////////////////////////////////////////////////
    //----------------------------------------------------
    if(olmTexture==1){
      var year = CY; //  
      if(CY <= 2018){
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        // Soil Texture - b0 = surface
        var olmTextureB0 = olm_texture.select('b0');
        // Soil Texture - b10 = 10cm
        var olmTextureB10 = olm_texture.select('b10');
        // Soil Texture - b30 = 30cm
        var olmTextureB30 = olm_texture.select('b30');
        // Soil Texture - b60 = 60cm
        var olmTextureB60 = olm_texture.select('b60');
        // Soil Texture - b100 = 100cm
        var olmTextureB100 = olm_texture.select('b100');
        // Soil Texture - b200 = 200cm
        var olmTextureB200 = olm_texture.select('b200');
        // export to google drive
        Export.image.toDrive({
          image: olmTextureB0,
          description: current_farm + '_olm_texture0cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmTextureB10,
          description: current_farm + '_olm_texture10cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmTextureB30,
          description: current_farm + '_olm_texture30cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmTextureB60,
          description: current_farm + '_olm_texture60cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmTextureB100,
          description: current_farm + '_olm_texture100cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmTextureB200,
          description: current_farm + '_olm_texture200cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // END IF YEAR <= 2018
    } // END TEXTURE
    
    ///////////////////////////////////////////////////////////
    // OPEN LAND MAP - BULK DENSITY - (250m)
    ///////////////////////////////////////////////////////////
    //----------------------------------------------------
    if(olmBulkDens==1){
      var year = CY; //  
      if(CY <= 2018){
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        // BULK DENSITY - b0 = surface
        var olmBulkdensityB0 = olm_bulkdensity.select('b0');
        // BULK DENSITY - b10 = 10cm
        var olmBulkdensityB10 = olm_bulkdensity.select('b10');
        // BULK DENSITY - b30 = 30cm
        var olmBulkdensityB30 = olm_bulkdensity.select('b30');
        // BULK DENSITY - b60 = 60cm
        var olmBulkdensityB60 = olm_bulkdensity.select('b60');
        // BULK DENSITY - b100 = 100cm
        var olmBulkdensityB100 = olm_bulkdensity.select('b100');
        // BULK DENSITY - b200 = 200cm
        var olmBulkdensityB200 = olm_bulkdensity.select('b200');
        // export to google drive
        Export.image.toDrive({
          image: olmBulkdensityB0,
          description: current_farm + '_olm_bulkdensity0cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmBulkdensityB10,
          description: current_farm + '_olm_bulkdensity10cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmBulkdensityB30,
          description: current_farm + '_olm_bulkdensity30cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmBulkdensityB60,
          description: current_farm + '_olm_bulkdensity60cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmBulkdensityB100,
          description: current_farm + '_olm_bulkdensity100cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmBulkdensityB200,
          description: current_farm + '_olm_bulkdensity200cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // END IF YEAR <= 2018
    } // END BULK DENSITY
    
    ///////////////////////////////////////////////////////////
    // OPEN LAND MAP - CLAY CONTENT - (250m)
    ///////////////////////////////////////////////////////////
    //----------------------------------------------------
    if(olmClayContent==1){
      var year = CY; //  
      if(CY <= 2018){
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        // CLAY CONTENT - b0 = surface
        var olmClayContentB0 = olm_claycontent.select('b0');
        // CLAY CONTENT - b10 = 10cm
        var olmClayContentB10 = olm_claycontent.select('b10');
        // CLAY CONTENT - b30 = 30cm
        var olmClayContentB30 = olm_claycontent.select('b30');
        // CLAY CONTENT - b60 = 60cm
        var olmClayContentB60 = olm_claycontent.select('b60');
        // CLAY CONTENT - b100 = 100cm
        var olmClayContentB100 = olm_claycontent.select('b100');
        // CLAY CONTENT - b200 = 200cm
        var olmClayContentB200 = olm_claycontent.select('b200');
        // export to google drive
        Export.image.toDrive({
          image: olmClayContentB0,
          description: current_farm + '_olm_claycontent0cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmClayContentB10,
          description: current_farm + '_olm_claycontent10cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmClayContentB30,
          description: current_farm + '_olm_claycontent30cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmClayContentB60,
          description: current_farm + '_olm_claycontent60cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmClayContentB100,
          description: current_farm + '_olm_claycontent100cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmClayContentB200,
          description: current_farm + '_olm_claycontent200cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // END IF YEAR <= 2018
    } // END CLAY CONTENT
    
    ///////////////////////////////////////////////////////////
    // OPEN LAND MAP - SAND CONTENT - (250m)
    ///////////////////////////////////////////////////////////
    //----------------------------------------------------
    if(olmSandContent==1){
      var year = CY; //  
      if(CY <= 2018){
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        // SAND CONTENT - b0 = surface
        var olmSandContentB0 = olm_sandcontent.select('b0');
        // SAND CONTENT - b10 = 10cm
        var olmSandContentB10 = olm_sandcontent.select('b10');
        // SAND CONTENT - b30 = 30cm
        var olmSandContentB30 = olm_sandcontent.select('b30');
        // SAND CONTENT - b60 = 60cm
        var olmSandContentB60 = olm_sandcontent.select('b60');
        // SAND CONTENT - b100 = 100cm
        var olmSandContentB100 = olm_sandcontent.select('b100');
        // SAND CONTENT - b200 = 200cm
        var olmSandContentB200 = olm_sandcontent.select('b200');
        // export to google drive
        Export.image.toDrive({
          image: olmSandContentB0,
          description: current_farm + '_olm_sandcontent0cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmSandContentB10,
          description: current_farm + '_olm_sandcontent10cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmSandContentB30,
          description: current_farm + '_olm_sandcontent30cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmSandContentB60,
          description: current_farm + '_olm_sandcontent60cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmSandContentB100,
          description: current_farm + '_olm_sandcontent100cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmSandContentB200,
          description: current_farm + '_olm_sandcontent200cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // END IF YEAR <= 2018
    } // END SAND CONTENT
    
    ///////////////////////////////////////////////////////////
    // OPEN LAND MAP - pH in WATER - (250m)
    ///////////////////////////////////////////////////////////
    //----------------------------------------------------
    if(olmpHw==1){
      var year = CY; //  
      if(CY <= 2018){
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        // pH(w) - b0 = surface
        var olmpHwB0 = olm_phw.select('b0');
        // pH(w) - b10 = 10cm
        var olmpHwB10 = olm_phw.select('b10');
        // pH(w) - b30 = 30cm
        var olmpHwB30 = olm_phw.select('b30');
        // pH(w) - b60 = 60cm
        var olmpHwB60 = olm_phw.select('b60');
        // pH(w) - b100 = 100cm
        var olmpHwB100 = olm_phw.select('b100');
        // pH(w) - b200 = 200cm
        var olmpHwB200 = olm_phw.select('b200');
        // export to google drive
        Export.image.toDrive({
          image: olmpHwB0,
          description: current_farm + '_olm_phw0cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmpHwB10,
          description: current_farm + '_olm_phw10cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmpHwB30,
          description: current_farm + '_olm_phw30cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmpHwB60,
          description: current_farm + '_olm_phw60cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmpHwB100,
          description: current_farm + '_olm_phw100cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmpHwB200,
          description: current_farm + '_olm_phw200cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // END IF YEAR <= 2018
    } // END pH(w) 
    
    ///////////////////////////////////////////////////////////
    // OPEN LAND MAP - WATER CONTENT @ 33kpa (Field Capacity) - (250m)
    ///////////////////////////////////////////////////////////
    //----------------------------------------------------
    if(olmWaterContent==1){
      var year = CY; //  
      if(CY <= 2018){
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        // WATER CONTENT - b0 = surface
        var olmWaterContentB0 = olm_watercontent.select('b0');
        // WATER CONTENT - b10 = 10cm
        var olmWaterContentB10 = olm_watercontent.select('b10');
        // WATER CONTENT - b30 = 30cm
        var olmWaterContentB30 = olm_watercontent.select('b30');
        // WATER CONTENT - b60 = 60cm
        var olmWaterContentB60 = olm_watercontent.select('b60');
        // WATER CONTENT - b100 = 100cm
        var olmWaterContentB100 = olm_watercontent.select('b100');
        // WATER CONTENT - b200 = 200cm
        var olmWaterContentB200 = olm_watercontent.select('b200');
        // export to google drive
        Export.image.toDrive({
          image: olmWaterContentB0,
          description: current_farm + '_olm_watercontent0cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmWaterContentB10,
          description: current_farm + '_olm_watercontent10cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmWaterContentB30,
          description: current_farm + '_olm_watercontent30cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmWaterContentB60,
          description: current_farm + '_olm_watercontent60cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmWaterContentB100,
          description: current_farm + '_olm_watercontent100cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmWaterContentB200,
          description: current_farm + '_olm_watercontent200cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // END IF YEAR <= 2018
    } // END WATER CONTENT
    
    ///////////////////////////////////////////////////////////
    // OPEN LAND MAP - ORGANIC CARBON CONTENT - (250m)
    ///////////////////////////////////////////////////////////
    //----------------------------------------------------
    if(olmCarbonContent==1){
      var year = CY; //  
      if(CY <= 2018){
        var mask = ee.FeatureCollection(farm_name);
        var mask_geometry = mask.geometry().bounds();
        // CARBON CONTENT - b0 = surface
        var olmCarbonContentB0 = olm_carboncontent.select('b0');
        // CARBON CONTENT - b10 = 10cm
        var olmCarbonContentB10 = olm_carboncontent.select('b10');
        // CARBON CONTENT - b30 = 30cm
        var olmCarbonContentB30 = olm_carboncontent.select('b30');
        // CARBON CONTENT - b60 = 60cm
        var olmCarbonContentB60 = olm_carboncontent.select('b60');
        // CARBON CONTENT - b100 = 100cm
        var olmCarbonContentB100 = olm_carboncontent.select('b100');
        // CARBON CONTENT - b200 = 200cm
        var olmCarbonContentB200 = olm_carboncontent.select('b200');
        // export to google drive
        Export.image.toDrive({
          image: olmCarbonContentB0,
          description: current_farm + '_olm_carboncontent0cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmCarbonContentB10,
          description: current_farm + '_olm_carboncontent10cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmCarbonContentB30,
          description: current_farm + '_olm_carboncontent30cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmCarbonContentB60,
          description: current_farm + '_olm_carboncontent60cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmCarbonContentB100,
          description: current_farm + '_olm_carboncontent100cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
        Export.image.toDrive({
          image: olmCarbonContentB200,
          description: current_farm + '_olm_carboncontent200cm_250m'  + '_' + year,
          folder:'GEE_Surfaces',
          scale: 250,
          crs: 'EPSG:4326',
          region: mask_geometry
        });
      } // END IF YEAR <= 2018
    } // END ORGANIC CARBON CONTENT
    //----------------------------------------------------
    
    //********************************************************************
    // END FARMS LOOP
    //********************************************************************
  }
  //********************************************************************
  // END YEARS LOOP
  //********************************************************************
}







