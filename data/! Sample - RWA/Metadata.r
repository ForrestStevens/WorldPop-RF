#####
##	BEGIN:	Standard configuration options

##	Configure the country abbreviation and name:
country_name <- "Rwanda"


##	These should be set before the Metadata.r script is run, included from
##		the processing scripts in the RF /src folder.  But if they aren't and
##		Metadata.r is run separately they should be configured here otherwise
##		the following configuration information will not run:
#country <- "RWA"
#rf_version = "2b"
#root_path <- "D:\\Documents\\Graduate School\\Research\\Population\\Analysis\\RF\\Working RF\\"
#source("01.0 - Configuration.py.r")


##	Configuration options for RandomForest modeling and prediction:

##  If we are using a set of covariates from another country set the 
##		fixed_set variable to specify which will then be used to fix the 
##		covariate list to an existing randomForest object, otherwise, 
##		use the full set from the covariate metadata if it is NULL (or set
##		to the "country" variable that you specify above).
##
##    Note that the fixed_set flag also changes the way that the 
##		randomForest object is created below by eliminating the variable 
##		elimination routine:
#fixed_set <- "VNM"
#fixed_set <- c("VNM", "KHM")
fixed_set <- NULL


##	Configure project and default data folders (typically these are static
##		and should never change as they are constructed here mainly to 
##		take into account possibly having moved the RF base folder):

##	Setup output and project paths, creating empty folders if necessary:
output_path <- paste(root_path, "/output/", country, "/", sep="")
dir.create(output_path, showWarnings=FALSE)
output_path_tmp <- paste(root_path, "/output/", country, "/tmp/", sep="")
dir.create(output_path_tmp, showWarnings=FALSE)

project_path <- paste(root_path, "/data/", country, "/", sep="")

##	END:	Standard configuration options
#####



#####
##	BEGIN:	Metadata descriptions:

##	During documentation these list items will be added describing the primary classes of datsets we use and processes that determine which covariates will be added to the modeling output:
metadata <- list()

##	Required datasets that should be placed and described by hand:
##		
##		Census, Landcover, NPP (created from "Data Preparation, R")
##
##	Default datasets that will be created from existing data sources if no
##		shapefile or image is found in the folder of the same name:
##		
##		Lights, Roads, Rivers, Waterbodies, Populated, Protected, Elevation, 
##		Temp, Precip
##
##	Any additional ancillary datasets can be specified as long as they have
##		a one word name and that name starts with three unique letters.  The
##		name is assumed to be the folder name but the file name of the 
##		shapefile or image can be anything you want it to be.


##	Required datasets, should be updated for every country either here or
##		overwritten in the custom ancillary datasets below:


var_name <- "Census"
var_folder <- paste(project_path, var_name, "/", sep="")
metadata[[var_name]] <- list(
	dataset_folder = var_name,
	dataset_title = "",
	dataset_source = "",
	dataset_name = list.files(var_folder, "shp$"),
	dataset_description = "Required fields for map production are ADMINID and ADMINPOP.",
	dataset_class = "polygon",
	derived = c("area", "buff", "zones")
)
metadata[[var_name]][["path"]] <- paste(var_folder, metadata[[var_name]]$dataset_name, sep="")


var_name <- "Landcover"
var_folder <- paste(project_path, var_name, "/", sep="")
metadata[[var_name]] <- list(
	dataset_folder = var_name,
	dataset_title = "",
	dataset_source = "",
	dataset_name = list.files(var_folder, "tif$|img$"),
	dataset_description = "",
	dataset_class = "raster",
	derived = c( 
		"prp011", "cls011", "dst011",
		"prp040", "cls040", "dst040",
		"prp130", "cls130", "dst130",
		"prp140", "cls140", "dst140",
		"prp150", "cls150", "dst150",
		"prp160", "cls160", "dst160",
		"prp190", "cls190", "dst190",
		"prp200", "cls200", "dst200",
		"prp210", "cls210", "dst210",
		"prp230", "cls230", "dst230",
		"prp240", "cls240", "dst240",
		"prp250", "cls250", "dst250",
		"prpBLT", "clsBLT", "dstBLT"
	)
)
metadata[[var_name]][["path"]] <- paste(var_folder, metadata[[var_name]]$dataset_name, sep="")



##	Default data that are generated from default data sources unless 
##		replaced by hand in the custom ancillary metadata below:


var_name <- "NPP"
metadata[[var_name]] <- list(
	dataset_folder = var_name,
	dataset_title = "MODIS 17A3 2010 Estimated Net Primary Productivity, 1km",
	dataset_source = "United States Geological Survey (USGS)",
	dataset_name = "DEFAULT: MODIS 17A3 2010",
	dataset_description = "MODIS 17A3 version-55 derived estimates of net primary productivity for the year 2010, estimated for 1km pixel sizes and subset and resampled to match the available land cover and final population map output requirements.",
	dataset_class = "raster",
	derived = c("")
)
metadata[[var_name]][["path"]] <- paste(project_path, metadata[[var_name]]$dataset_folder, "/", metadata[[var_name]]$dataset_name, sep="")


var_name <- "Lights"
metadata[[var_name]] <- list(
	dataset_folder = var_name,
	dataset_title = "Suomi NPP VIIRS-Derived 2012 Lights at Night, 15 arc-second",
	dataset_source = "http://ngdc.noaa.gov/eog/viirs/download_viirs_ntl.html",
	dataset_name = "DEFAULT: VIIRS 2012",
	dataset_description = "These 'Lights at Night' data were derived from imagery collected by the Suomi National Polar-orbiting Partnership (NPP) Visible Infrared Imaging Radiometer Suite (VIIRS) sensor.  Data were collected in 2012 on moonless nights and though background noise associated with fires, gas-flares, volcanoes or aurora have not been removed it represents the best-available data for night-time light production.",
	dataset_class = "raster",
	derived = c("")
)
metadata[[var_name]][["path"]] <- paste(project_path, metadata[[var_name]]$dataset_folder, "/", metadata[[var_name]]$dataset_name, sep="")


var_name <- "Temp"
metadata[[var_name]] <- list(
	dataset_folder = var_name,
	dataset_title = "WorldClim/BioClim Mean Annual Temperature 1950-2000, 30 arc-second",
	dataset_source = "http://www.worldclim.org/current",
	dataset_name = "DEFAULT: BIO1",
	dataset_description = "WorldClim/BioClim 1950-2000 mean annual precipitation (BIO12) and mean annual temperature (BIO1) estimates (Hijmans et al., 2005) were downloaded, mosaicked and subset to match the extent of our land cover data for the mapping of this region.",
	dataset_class = "raster",
	derived = c("")
)
metadata[[var_name]][["path"]] <- paste(project_path, metadata[[var_name]]$dataset_folder, "/", metadata[[var_name]]$dataset_name, sep="")


var_name <- "Precip"
metadata[[var_name]] <- list(
	dataset_folder = var_name,
	dataset_title = "WorldClim/BioClim Mean Annual Precipitation 1950-2000, 30 arc-second",
	dataset_source = "http://www.worldclim.org/current",
	dataset_name = "DEFAULT: BIO12",
	dataset_description = "WorldClim/BioClim 1950-2000 mean annual precipitation (BIO12) and mean annual temperature (BIO1) estimates (Hijmans et al., 2005) were downloaded, mosaicked and subset to match the extent of our land cover data for the mapping of this region.",
	dataset_class = "raster",
	derived = c("")
)
metadata[[var_name]][["path"]] <- paste(project_path, metadata[[var_name]]$dataset_folder, "/", metadata[[var_name]]$dataset_name, sep="")


var_name <- "Roads"
metadata[[var_name]] <- list(
	dataset_folder = var_name,
	dataset_title = "Road Network",
	dataset_source = "National Geospatial-Intelligence Agency (NGA), http://geoengine.nga.mil/geospatial/SW_TOOLS/NIMAMUSE/webinter/rast_roam.html",
	dataset_name = "DEFAULT: trans/roadl",
	dataset_description = "The VMAP0 data area downloaded as separate files, grouped roughly by continent, and merged into individual shapefiles for subsetting and further processing for population mapping efforts.  These data were obtained directly from the original VMAP0 data sources provided by the NGA and pre-processed using Military Analyst in ArcGIS 10.0.",
	dataset_class = "linear",
	derived = c("dst")
)
metadata[[var_name]][["path"]] <- paste(project_path, metadata[[var_name]]$dataset_folder, "/", metadata[[var_name]]$dataset_name, sep="")


var_name <- "Rivers"
metadata[[var_name]] <- list(
	dataset_folder = var_name,
	dataset_title = "River Network",
	dataset_source = "National Geospatial-Intelligence Agency (NGA), http://geoengine.nga.mil/geospatial/SW_TOOLS/NIMAMUSE/webinter/rast_roam.html",
	dataset_name = "DEFAULT: hydro/watrcrsl",
	dataset_description = "The VMAP0 data area downloaded as separate files, grouped roughly by continent, and merged into individual shapefiles for subsetting and further processing for population mapping efforts.  These data were obtained directly from the original VMAP0 data sources provided by the NGA and pre-processed using Military Analyst in ArcGIS 10.0.",
	dataset_class = "linear",
	derived = c("dst")
)
metadata[[var_name]][["path"]] <- paste(project_path, metadata[[var_name]]$dataset_folder, "/", metadata[[var_name]]$dataset_name, sep="")


var_name <- "Populated"
metadata[[var_name]] <- list(
	dataset_folder = var_name,
	dataset_title = "Populated Places",
	dataset_source = "National Geospatial-Intelligence Agency (NGA), http://geoengine.nga.mil/geospatial/SW_TOOLS/NIMAMUSE/webinter/rast_roam.html",
	dataset_name = "DEFAULT: Merged pop/builtupp, pop/builtupa, pop/mispopp", 
	dataset_description = "The VMAP0 data area downloaded as separate files, grouped roughly by continent, and merged into individual shapefiles for subsetting and further processing for population mapping efforts.  These data were obtained directly from the original VMAP0 data sources provided by the NGA and pre-processed using Military Analyst in ArcGIS 10.0.  Point data sources are buffered to 100 m and then all polygon data sources are merged to a single shapefile prior to processing.",
	dataset_class = "polygon",
	derived = c("cls", "dst", "prp")
)
metadata[[var_name]][["path"]] <- paste(project_path, metadata[[var_name]]$dataset_folder, "/", metadata[[var_name]]$dataset_name, sep="")


var_name <- "Waterbodies"
metadata[[var_name]] <- list(
	dataset_folder = var_name,
	dataset_title = "Inland Waterbodies",
	dataset_source = "National Geospatial-Intelligence Agency (NGA), http://geoengine.nga.mil/geospatial/SW_TOOLS/NIMAMUSE/webinter/rast_roam.html",
	dataset_name = "DEFAULT: hydro/watrcrsl",
	dataset_description = "The VMAP0 data area downloaded as separate files, grouped roughly by continent, and merged into individual shapefiles for subsetting and further processing for population mapping efforts.  These data were obtained directly from the original VMAP0 data sources provided by the NGA and pre-processed using Military Analyst in ArcGIS 10.0.",
	dataset_class = "polygon",
	derived = c("cls", "dst", "prp")
)
metadata[[var_name]][["path"]] <- paste(project_path, metadata[[var_name]]$dataset_folder, "/", metadata[[var_name]]$dataset_name, sep="")


var_name <- "Protected"
metadata[[var_name]] <- list(
	dataset_folder = var_name,
	dataset_title = "Protected Areas",
	dataset_source = "World Database on Protected Areas, Downloaded September, 2012, UNEP, http://www.wdpa.org, http://protectedplanet.net",
	dataset_name = "DEFAULT: WDPAfgdb_Sept2012.gdb",
	dataset_description = "These data are compiled by UNEP and distributed via the Protected Planet website.  All protected areas were downloaded regardless of International Union for Conservation of Nature (IUCN) or any other designation, so they include sanctuaries, national parks, game reserves, World Heritage Sites, etc.",
	dataset_class = "polygon",
	derived = c("cls", "dst", "prp")
)
metadata[[var_name]][["path"]] <- paste(project_path, metadata[[var_name]]$dataset_folder, "/", metadata[[var_name]]$dataset_name, sep="")


var_name <- "Urban"
metadata[[var_name]] <- list(
	dataset_folder = var_name,
	dataset_title = "Urban Extents",
	dataset_source = "Schneider, et al., United Nations",
	dataset_name = "DEFAULT: schneider-urban.shp",
	dataset_description = "These data were constructed from MODIS-derived imagery and provided to WorldPop researchers by Schneider, et al. as part of a global urban extents datasets.",
	dataset_class = "polygon",
	derived = c("cls", "dst", "prp")
)
metadata[[var_name]][["path"]] <- paste(project_path, metadata[[var_name]]$dataset_folder, "/", metadata[[var_name]]$dataset_name, sep="")


var_name <- "Elevation"
metadata[[var_name]] <- list(
	dataset_folder = var_name,
	dataset_title = "Elevation and Derived Slope, 3 second",
	dataset_source = "HydroSHEDS Void-Filled DEM (Lehnert, et al., 2006), http://hydrosheds.cr.usgs.gov/dataavail.php",
	dataset_name = "DEFAULT: Void-Filled DEM.gdb",
	dataset_description = "The HydroSHEDS data are the result of an effort to provide a globally consistent dataset consisting of NASA's Shuttle Radar Topography Mission (SRTM) data and have been processed, void-filled and corrected for use at large scales.",
	dataset_class = "raster",
	derived = c("", "slope")
)
metadata[[var_name]][["path"]] <- paste(project_path, metadata[[var_name]]$dataset_folder, "/", metadata[[var_name]]$dataset_name, sep="")



##	All "Custom Ancillary" datasets should appear below.  They may replace
##		those named above but make sure that if they are not replacements
##		that their var_name is a single word and begins with three letters
##		that have not already been used:


#var_name <- "Roads"
#var_folder <- paste(project_path, var_name, "/", sep="")
#metadata[[var_name]] <- list(
#	dataset_folder = var_name,
#	dataset_title = "Roads (UN WFP)",
#	dataset_source = "UN World Food Programme, 2013, http://www.wfp.org",
#	dataset_name = list.files(var_folder, "shp$"),
#	dataset_description = "These data were downloaded as part of a per-country package of data layers made availalble as shapefiles through the World Food Programme website.",
#	dataset_class = "linear",
#	derived = c("dst")
#)
#metadata[[var_name]][["path"]] <- paste(project_path, metadata[[var_name]]$dataset_folder, "/", metadata[[var_name]]$dataset_name, sep="")


##	Alternatively we can specify them in a vector and build the metadata 
##		structures programmatically and optionally overwriting any of the 
##		default data sources above:


var_names <- c(
	"Census",
	"Landcover",
	"Buildings",
  "Points",
  "Rivers",
  "Roads",
  "Uses"
)

dataset_classes <- c(
	"polygon",
	"raster",
	"polygon",
  "point",
	"linear",
	"linear",
	"polygon"
)

dataset_titles <- c(
	"Rwandan Census, 2002",
	"Remotely-sensed, Classified Landcover",
	"Digitized Building Locations (OSM), 2013",
	"Points of Interest Locations (OSM), 2013",
	"Rivers (OSM), 2013",
	"Roads (OSM), 2013",
	"Delineated Land Uses (OSM), 2013"
)

dataset_sources <- c(
	"National Institute of Statistics of Rwanda, 2002",
	"Fused Globcover 2009 (300m) With Landsat-Derived Urban/Rural Cover (100m)",
	"Open Street Map, Downloaded 2013-09-16, http://extract.bbbike.org/",
	"Open Street Map, Downloaded 2013-09-16, http://extract.bbbike.org/",
	"Open Street Map, Downloaded 2013-09-16, http://extract.bbbike.org/",
	"Open Street Map, Downloaded 2013-09-16, http://extract.bbbike.org/",
	"Open Street Map, Downloaded 2013-09-16, http://extract.bbbike.org/"
)

dataset_descriptions <- c(
	"These high spatial resolution census block data were attained through in-country partners for 2002.",
	"Land cover information was combined from a GlobCover 2009 coverage and fused with Landsat-derived urban/rural built area classification to construct a single land cover dataset.",
	"These data were downloaded as part of a per-country package of data layers made availalble as shapefiles through the http://extract.bbbike.org website, extracted from the Open Street Map (OSM) database.",
	"These data were downloaded as part of a per-country package of data layers made availalble as shapefiles through the http://extract.bbbike.org website, extracted from the Open Street Map (OSM) database.",
	"These data were downloaded as part of a per-country package of data layers made availalble as shapefiles through the http://extract.bbbike.org website, extracted from the Open Street Map (OSM) database.",
	"These data were downloaded as part of a per-country package of data layers made availalble as shapefiles through the http://extract.bbbike.org website, extracted from the Open Street Map (OSM) database.",
	"These data were downloaded as part of a per-country package of data layers made availalble as shapefiles through the http://extract.bbbike.org website, extracted from the Open Street Map (OSM) database."
)

dataset_derived <- list(
	c("area", "buff", "zones"),
	c( 
		"prp011", "cls011", "dst011",
		"prp040", "cls040", "dst040",
		"prp130", "cls130", "dst130",
		"prp140", "cls140", "dst140",
		"prp150", "cls150", "dst150",
		"prp160", "cls160", "dst160",
		"prp190", "cls190", "dst190",
		"prp200", "cls200", "dst200",
		"prp210", "cls210", "dst210",
		"prp230", "cls230", "dst230",
		"prp240", "cls240", "dst240",
		"prp250", "cls250", "dst250",
		"prpBLT", "clsBLT", "dstBLT"
	),
	c("prp", "cls", "dst"),
	c("prp", "cls", "dst"),
	c("prp", "cls", "dst"),
	c("prp", "cls", "dst"),
	c("prp", "cls", "dst")
)

for (i in 1:length(var_names)) {
	var_name <- var_names[i]
	var_folder <- paste(project_path, var_name, "/", sep="")

	metadata[[var_name]] <- list(
		dataset_folder = var_name,
		dataset_title = dataset_titles[i],
		dataset_source = dataset_sources[i],
		dataset_name = list.files(var_folder, "shp$|tif$|img$"),
		dataset_description = dataset_descriptions[i],
		dataset_class = dataset_classes[i],
		derived = dataset_derived[[i]]
	)
	metadata[[var_name]][["path"]] <- paste(var_folder, metadata[[var_name]]$dataset_name, sep="")
}


##	END:	Metadata descriptions:
#####
