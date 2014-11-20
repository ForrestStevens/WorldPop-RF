#####
##	BEGIN:	Set per-country configuration options

##	Parse main configuration file, which will set the country and root_path
##		variables:
execfile("01.0 - Configuration.py.r")


##	Output projection is configured for each country based on the most
##		appropriate output for distance/area considerations.  To find the
##		string appropriate for this string you should open up the appropriate
##		projection file (*.prj) in the:
##				C:\Program Files (x86)\ArcGIS\Desktop10.0\Coordinate Systems\
##		folder and paste its contents into a string below (or hand edit it
##		in case a UTM zone needs to be edited, etc.:

##	For KHM, VNM:
#intermediate_prj = "PROJCS['WGS_1984_UTM_Zone_48N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',105.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]"
##	For KEN:
#intermediate_prj = "PROJCS['WGS_1984_UTM_Zone_37N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',39.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]"
##	For NGA:
#intermediate_prj = "PROJCS['WGS_1984_UTM_Zone_32N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',9.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0],AUTHORITY['EPSG',32632]]"
##	For RWA, right on edge of 35S and 36S so modified to 35.5S with meridian 30.0:
intermediate_prj = "PROJCS['WGS_1984_UTM_Zone_35.5S',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',10000000.0],PARAMETER['Central_Meridian',30.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0],AUTHORITY['EPSG',32735.5]]"
##	For IDN:
#intermediate_prj = 'PROJCS["Asia_Lambert_Conformal_Conic",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",105.0],PARAMETER["Standard_Parallel_1",30.0],PARAMETER["Standard_Parallel_2",62.0],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0],AUTHORITY["ESRI",102012]]'
##	For PAN:
#intermediate_prj = 'PROJCS["Panama-Colon 1911 / Panama Lambert",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",294865.303],PARAMETER["Central_Meridian",-80.0],PARAMETER["Standard_Parallel_1",8.25],PARAMETER["Standard_Parallel_2",8.25],PARAMETER["Scale_Factor",0.99989909],PARAMETER["Latitude_Of_Origin",8.25],UNIT["Meter",1.0]]'
##	For CRI:
#intermediate_prj = 'PROJCS["WGS_1984_UTM_Zone_16.5N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-84.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0],AUTHORITY["EPSG",32616]]'
##	For NIC:
#intermediate_prj = 'PROJCS["WGS_1984_UTM_Zone_16.25N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-85.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0],AUTHORITY["EPSG",32616]]'
##	For SUR:
#intermediate_prj = "PROJCS['Zanderij_1972_UTM_Zone_21N',GEOGCS['GCS_Zanderij',DATUM['D_Zanderij',SPHEROID['International_1924',6378388.0,297.0]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-57.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0],AUTHORITY['EPSG',31121]]"
#intermediate_prj = "PROJCS['WGS_1984_UTM_Zone_21N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-57.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0],AUTHORITY['EPSG',32621]]"
##	For HAI:
#intermediate_prj = "PROJCS['WGS_1984_UTM_Zone_18.5N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-72.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0],AUTHORITY['EPSG',32618]]"

##	Alternative census folder:
##		This is a census folder other than the one used for the RF model
##		and predict_density.img file.  If the "! New Census" folder exists
##		then the census	file in this folder overrides the one in the /Census
##		folder during final map creation for census counts.  This folder must
##		be set to "! New Census/" because the reporting procedures need to
##		be able to tell if a newer census was used or not.  Its derived
##		products only consist of a reprojected dataset.  This is useful
##		if you'll using a much more detailed, but older census file to
##		generate RF model	predictions and want to use more concurrent data
##		for population	estimates.
##	NOTE: Setting or changing the below folder name will not affect anything
##		about how this folder is processed.  It will be processed if it exists
##		and should not be used if you don't intend for this behavior!!!
##	TODO:  Eventually this needs to be pulled from the Metadata.R file
##		instead of a hard coded directory presence.  If this directory exists
##		then the rest of the scripts (map production and metadata report
##		will assume it needs to be used.
#census_folder = "! New Census"


##	Should we skip processing and creation for any existing data sets:
skip_existing = True


##	END:	Set per-country configuration options
#####



#####
##	NOTICE:  In practice nothing below this line should need to be regularly
##		edited unless very specific details need to be changed about the
##		modeling process (e.g. speciyfing an extent or something detailed
##		about the Python processing).
#####



#####
##	BEGIN:	Set general configuration options

data_path = root_path + "/data/"
output_path = root_path + "/output/"

VMAP0_path = data_path + "/! Defaults/VMAP0/"
WorldClim_path = data_path + "/! Defaults/WorldClim/"
HydroSHEDS_path = data_path + "/! Defaults/HydroSHEDS/"
WDPA_path = data_path + "/! Defaults/Protected Areas/"
UrbExtents_path = data_path + "/! Defaults/Urban Extents/"

#NightLights_path = data_path + "/! Defaults/NASA Night Lights/"
NightLights_path = data_path + "/! Defaults/NASA Night Lights, Unprocessed/"

##	END:	Set general configuration options
#####



#####
##	BEGIN: Import packages and set GeoProcessing environment

##	Import the Spatial Analyst extension for slope and raster algebra
##		functionality:
import os
import glob
import json
#from pprint import pprint

import arcpy
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *


##	Configure arcpy:

# Set the compression environment to LZ77:
arcpy.env.compression = "LZ77"

##	Should we overwrite any already existing derived data:
overwrite = True
arcpy.env.overwriteOutput = overwrite

##	Set our default to not build pyramids:
arcpy.env.pyramid = "PYRAMIDS 0"


##	END: Import packages and set GeoProcessing environment
#####



#####
##	BEGIN: Define utility functions

##	Create a directory function to check for and create data directories if they don't exist:
def ensure_dir(d):
	if not os.path.isdir(d):
		os.makedirs(d)
	return d


def process_linear(dataset_folder):
	##	Create distance rasters for each of our feature classes of interest:

	dataset_folder = dataset_folder
	dataset_name = dataset_folder.lower()

	print("PROCESSING:  " + dataset_name + " Distance")
	out_path = data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + "_dst" +".tif"

	if not os.path.isfile(out_path) or not skip_existing:
		outRas = arcpy.sa.EucDistance(data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + ".shp")
		#outRas.save(out_path)

		##	NOTE: Instead of just saving the output to a file using the arcpy 
		##		raster object method, we use CopyRaster here instead which allows 
		##		data compression to be implemented correctly.  We also set the 
		##		arcpy.env.pyramid variable to not produce pyramids but by default
		##		the CopyRaster call creates them, so we have to call the function
		##		to remove them according to the environment setting above:
		outRas = arcpy.CopyRaster_management(outRas, out_path, "DEFAULTS", "0", "-999", "NONE", "NONE", "32_BIT_FLOAT")
		outRas = arcpy.BuildPyramids_management(outRas)

		##	NOTE: So there's a bug in the arcpy raster optimization that
		##		often, even though you're specifying the .save() option will not
		##		actually write it out to disk until it thinks that you're really
		##		removing it from memory.  To force this to happen across versions
		##		of ArcGIS and arcpy set the raster object to None like this.
		##		Discovered here:
		##			http://gis.stackexchange.com/questions/46897/saving-rasters-in-a-python-for-loop-fails-only-on-last-iteration
		outRas = None


def process_point_area(dataset_folder):
	##	Process select feature data into binary masks, distance-to
	##		and proportions per	81 hectares (11 x 11 cell moving window,
	##		specified as a NbrCircle with radius 5):

	dataset_folder = dataset_folder
	dataset_name = dataset_folder.lower()

	print("PROCESSING:  " + dataset_folder + " Features")

	##	Use FID to create output raster from our features:
	##	NOTE:  This is one more kludge-y workaround for ArcGIS and arcpy
	##		bugginess.  For most cases except very large rasters like when
	##		dealing with China or Indonesia using one raster file for the
	##		post-feater-to-raster conversion seems to work.  But Python just
	##		bails for some large raster-cases.  Therefore we have to use this
	##		series of tmpRas objects in addition to the setting of things to
	##		None for the post-conversion process.  Hopefully one of these days
	##		we can remove this workaround...
	##			http://forums.arcgis.com/threads/81281-quot-python.exe-has-stopped-working-quot-on-raster.save
	tmpRas1 = data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + "_FID_tmp.tif"
	tmpRas2 = data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + "_FID.tif"
	tmpRas = data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + "_FID.tif"

	if not os.path.isfile(tmpRas2) or not skip_existing:
		##	See note below about PolygonToRaster and it not respecting the
		##		snap-to raster environment... Same here with all FeatureToRaster
		##		calls...

		##	Also, in ArcGIS 10.1 there's a bug with in FeatureRoRaster and
		##		other raster conversion algorithms that creates blocks of 0
		##		data in areas that occur outside the set processing extent.  I've
		##		found only one workaround and that's to remove the processing
		##		extent, perform the conversion, and then reset the processing
		##		extent, after which do a simple calculation on the converted
		##		raster and save that as our output.  It's quite gankty but for
		##		10.1 I don't know of another workaround.  So instead of:
		#arcpy.FeatureToRaster_conversion(data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + ".shp", "FID", tmpRas, 100)

		##		we do this:
		tmp_extent = arcpy.env.extent

		##	This line can be commented to turn this behavior off (for example if you
		##		end up with no data bars in ArcGIS 10.0):
		#arcpy.env.extent = ""

		outRas = arcpy.FeatureToRaster_conversion(data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + ".shp", "FID", tmpRas1, 100)
		print("PROCESSED:  " + dataset_folder + " Feature to Raster")
		arcpy.env.extent = tmp_extent
		outRas = None
		outRas = Raster(tmpRas1) * 1
		print("PROCESSED:  " + dataset_folder + " Multiplication")
		#outRas.save( tmpRas2 )
		outRas = arcpy.CopyRaster_management(outRas, tmpRas2, "DEFAULTS", "0", "4294967295", "NONE", "NONE", "32_BIT_UNSIGNED")
		outRas = arcpy.BuildPyramids_management(outRas)
		outRas = None
		print("PROCESSED:  " + dataset_folder + " Output Saved")


		##	Create distance-to raster:
		out_path = data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + "_dst.tif"
		outRas = arcpy.sa.EucDistance(tmpRas)
		#outRas.save(out_path)
		outRas = arcpy.CopyRaster_management(outRas, out_path, "DEFAULTS", "0", "-999", "NONE", "NONE", "32_BIT_FLOAT")
		outRas = arcpy.BuildPyramids_management(outRas)
		outRas = None

		##	Save the binary classification raster:
		out_path = data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + "_cls.tif"
		outRas =  IsNull(Raster(tmpRas)) == 0
		#outRas.save(out_path)
		outRas = arcpy.CopyRaster_management(outRas, out_path, "DEFAULTS", "0", "3", "NONE", "NONE", "2_BIT")
		outRas = arcpy.BuildPyramids_management(outRas)
		outRas = None

		###	Calculate proportion of cover:
		#outRas = Raster(out_path)
		#out_path = data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + "_prp.tif"
		##outRas = arcpy.sa.FocalStatistics(outRas,NbrRectangle(11,11,"CELL"),"MEAN","DATA")
		#outRas = arcpy.sa.FocalStatistics(outRas,NbrCircle(5,"CELL"),"MEAN","DATA")
		##outRas.save(out_path)
		#outRas = arcpy.CopyRaster_management(outRas, out_path, "DEFAULTS", "0", "-999", "NONE", "NONE", "32_BIT_FLOAT")
		#outRas = arcpy.BuildPyramids_management(outRas)
		#outRas = None


def process_raster_binary(dataset_folder):
	dataset_folder = dataset_folder
	dataset_name = dataset_folder.lower()

	#	Clip and project our built raster:
	print("PROCESSING:  " + dataset_folder)

	#	Instead of hardcoding it above we'll look for a custom TIF or IMG
	#		file that should be in the directory:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.img" ) + glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.tif" ))
	if dataset_path:
		input_name = os.path.basename( dataset_path[0] )
	else:
		input_name = ""

	in_path = data_path + country + "/" + dataset_folder + "/" + input_name

	output_name = dataset_name + "_cls.tif"

	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	data_desc = arcpy.Describe(in_path)
	input_prj = data_desc.SpatialReference.exportToString()


	if not os.path.isfile(out_path) or not skip_existing:
		arcpy.ProjectRaster_management(in_path,out_path,intermediate_prj,"NEAREST","100","#","#",input_prj)

		#print("Proportion:  " + dataset_folder)
		##outRas = arcpy.sa.FocalStatistics(out_path,NbrRectangle(11,11,"CELL"),"MEAN","DATA")
		#outRas = arcpy.sa.FocalStatistics(out_path,NbrCircle(5,"CELL"),"MEAN","DATA")

		#out_path = out_path[:-8] + "_prp.tif"
		##outRas.save(out_path)
		#outRas = arcpy.CopyRaster_management(outRas, out_path, "DEFAULTS", "0", "-999", "NONE", "NONE", "32_BIT_FLOAT")
		#outRas = arcpy.BuildPyramids_management(outRas)
		#outRas = None


		print("Distance:  " + dataset_folder)
		lcRas = Raster(out_path[:-8] + "_cls.tif")
		tmpRas =  SetNull(lcRas != 1, 1)
		outRas = arcpy.sa.EucDistance(tmpRas)

		out_path = out_path[:-8] + "_dst.tif"
		#outRas.save(out_path)
		outRas = arcpy.CopyRaster_management(outRas, out_path, "DEFAULTS", "0", "-999", "NONE", "NONE", "32_BIT_FLOAT")
		outRas = arcpy.BuildPyramids_management(outRas)
		outRas = None


def process_raster_continuous(dataset_folder):
	dataset_folder = dataset_folder
	dataset_name = dataset_folder.lower()

	#	Clip and project our built raster:
	print("PROCESSING:  " + dataset_folder)

	#	Instead of hardcoding it above we'll look for a custom TIF or IMG
	#		file that should be in the directory:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.img" ) + glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.tif" ))
	if dataset_path:
		input_name = os.path.basename( dataset_path[0] )
	else:
		print("ERROR:  No .tif or .img file was found in the " + dataset_folder + " data folder!")
		exit()

	in_path = data_path + country + "/" + dataset_folder + "/" + input_name

	output_name = dataset_name + ".tif"

	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	data_desc = arcpy.Describe(in_path)
	input_prj = data_desc.SpatialReference.exportToString()

	if not os.path.isfile(out_path) or not skip_existing:
		arcpy.ProjectRaster_management(in_path,out_path,intermediate_prj,"NEAREST","100","#","#",input_prj)


##	END: Define utility functions
#####



#####
##	BEGIN: Data pre-processing for default datasets

##	Ensure that we have output and temporary directory for this country:
tmp_path = ensure_dir(output_path)
tmp_path = ensure_dir(output_path + country)
tmp_path = ensure_dir(output_path + country + "/tmp")
tmp_path = output_path + country + "/tmp/"


###	Now that we have a tmp_path specified, load the JSON objects that
###		specify the metadata and covariates specific to each dataset we
###		need to process:
#metadata_json = open(tmp_path + "metadata.json")
#covariates_json = open(tmp_path + "covariates.json")
#
#metadata = json.load(metadata_json)
#covariates = json.load(covariates_json)


##	TODO: I need to configure the covariate processing to respect the
##		covariate and processing flags for each metadata section.  Right now
##		all default outputs are generated depending on the output type
##		(ShapeFile, binary raster, continuous raster, etc.):
#pprint(metadata)
#metadata["Census"]
#covariates["Census"]
#metadata["Census"]["dataset_name"]


##	Get a list of all existing folders in the country's data folder:
dataset_folders = os.walk(data_path + country + "/").next()[1]


##	Now process all default datasets or datasets that require some special
##		handling by first checking to see if the folder exists, assuming
##		that if it does it has been intentionally created by the user or
##		the "Data Preparation, R.r" script as specified in the Metadata.r
##		file:

##	First we check to see if the "Census" folder exists, and if not throw
##		an error and exit().  If it does, remove it from the list of
##		directories to process and proceed:
if not ("Census" in dataset_folders):
	print("ERROR:  No \"Census\" folder found!  This is required and indicates you possibly did not run the \"Data Preparation, R.r\" script or specify configuration options correctly in this Python processing script!")
	exit()
else :
	dataset_folders.remove("Census")

##	Project census data:
dataset_folder = "Census"
print("PROCESSING:  " + dataset_folder)

tmp_path = ensure_dir(data_path + country + "/" + dataset_folder)

##	Instead of hardcoding it above we'll just pull in the only shapefile
##		that should be in the directory:
dataset_path = glob.glob( tmp_path + "/" + "*.shp" )
if dataset_path:
	dataset_name = os.path.basename( dataset_path[0] )
else:
	print("ERROR:  No census file found!")
	exit()

in_path = tmp_path + "/" + dataset_name

output_name = "census.shp"
tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
out_path = tmp_path + "/" + output_name

data_desc = arcpy.Describe(in_path)
input_prj = data_desc.SpatialReference.exportToString()

if not os.path.isfile(out_path) or not skip_existing:
	arcpy.Project_management(in_path, out_path, intermediate_prj,"#",input_prj)


##	Calculate polygon area now that we're projected:
in_path = out_path

output_name = "census_area.shp"
out_path = data_path + country + "/" + dataset_folder + "/Derived/" + output_name

if not os.path.isfile(out_path) or not skip_existing:
	arcpy.CalculateAreas_stats(in_path,out_path)


##	Create buffer for processing:
in_path = out_path

output_name = "census_buffer.shp"
buffer_path = data_path + country + "/" + dataset_folder + "/Derived/" + output_name

if not os.path.isfile(buffer_path) or not skip_existing:
	arcpy.Buffer_analysis(in_path, buffer_path,"10 Kilometers","FULL","ROUND","ALL","#")



##	Set the output coordinate system and extent to that of our buffer:
arcpy.env.outputCoordinateSystem = buffer_path

##	NOTE: To properly set the extent we need to query the description of
##		the shapefile by hand, for some reason setting the extent directly
##		to the name of a dataset doesn't work...
desc = arcpy.Describe(buffer_path)
arcpy.env.extent = desc.Extent



##	Now we need to process an alternative census folder if specified
##		above.  This is used in the case of having a finely detailed
##		census file located in our /Census folder with which we want to
##		estimate our RF model, but we want to redistribute counts based
##		on a coarser census file:
if ("! New Census" in dataset_folders):
	dataset_folder = "! New Census"

	##	Project census data:
	print("PROCESSING:  " + dataset_folder)

	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder)

	##	Instead of hardcoding it above we'll just pull in the only shapefile
	##		that should be in the directory:
	dataset_path = glob.glob( tmp_path + "/" + "*.shp" )
	if dataset_path:
		dataset_name = os.path.basename( dataset_path[0] )
	else:
		print("ERROR:  No census file found!")
		exit()

	in_path = tmp_path + "/" + dataset_name

	output_name = "census.shp"
	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	data_desc = arcpy.Describe(in_path)
	input_prj = data_desc.SpatialReference.exportToString()

	if not os.path.isfile(out_path) or not skip_existing:
		tmp_extent = arcpy.env.extent
		arcpy.env.extent = ""
		arcpy.Project_management(in_path, out_path, intermediate_prj,"#",input_prj)
		arcpy.env.extent = tmp_extent



if not ("Landcover" in dataset_folders):
	print("ERROR:  No \"Landcover\" folder found!  This is required and indicates you possibly did not run the \"Data Preparation, R.r\" script or specify configuration options correctly in this Python processing script!")
	exit()
else :
	dataset_folders.remove("Landcover")

##	Clip and project our land cover raster:
dataset_folder = "Landcover"
print("PROCESSING:  " + dataset_folder)

##	Instead of hardcoding it above we'll just pull in the only TIF or IMG
##		file that should be in the directory:
dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.img" ) + glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.tif" ))
if dataset_path:
	dataset_name = os.path.basename( dataset_path[0] )
else:
	print("ERROR:  No land cover file found!")
	exit()

in_path = data_path + country + "/" + dataset_folder + "/" + dataset_name

output_name = "landcover.tif"

tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
landcover_path = tmp_path + "/" + output_name

data_desc = arcpy.Describe(in_path)
input_prj = data_desc.SpatialReference.exportToString()


if not os.path.isfile(landcover_path) or not skip_existing:
	arcpy.ProjectRaster_management(in_path,landcover_path,intermediate_prj,"NEAREST","100","#","#",input_prj)



##	Set our output extent and snap raster environment settings to match
##		that of our land cover:
arcpy.env.outputCoordinateSystem = landcover_path

desc = arcpy.Describe(landcover_path)
arcpy.env.extent = desc.Extent

arcpy.env.snapRaster = landcover_path



#####
##	 You can start commenting here if you need to subset the processing
##		in any way:
#####



##	Process landcover data into binary landcover masks and proportions per
##		81 hectares (11 x 11 cell moving window):
for lc in ['011', '040', '130', '140', '150', '160', '190', '200', '210', '230', '240', '250']:
	print(lc)
	out_path = landcover_path[:-4] + "_cls" + lc + ".tif"
	if not os.path.isfile(out_path) or not skip_existing:
		outRas =  Raster(landcover_path) == int(lc)
		#outRas.save(out_path)
		outRas = arcpy.CopyRaster_management(outRas, out_path, "DEFAULTS", "0", "3", "NONE", "NONE", "2_BIT")
		outRas = arcpy.BuildPyramids_management(outRas)
		outRas = None

	in_path = out_path
	#out_path = landcover_path[:-4] + "_prp" + lc + ".tif"
	#if not os.path.isfile(out_path) or not skip_existing:
	#	#outRas = arcpy.sa.FocalStatistics(in_path,NbrRectangle(11,11,"CELL"),"MEAN","DATA")
	#	outRas = arcpy.sa.FocalStatistics(in_path,NbrCircle(5,"CELL"),"MEAN","DATA")
	#	outRas.save(out_path)
	#	outRas = None


##	Process land cover combinations of interest:

##	Combined built class (land cover 190 and 240):
lc = "BLT"
print(lc)
out_path = landcover_path[:-4] + "_cls" + lc + ".tif"
if not os.path.isfile(out_path) or not skip_existing:
	outRas =  (Raster(landcover_path) == 190) + (Raster(landcover_path) == 240)
	#outRas.save(out_path)
	outRas = arcpy.CopyRaster_management(outRas, out_path, "DEFAULTS", "0", "3", "NONE", "NONE", "2_BIT")
	outRas = arcpy.BuildPyramids_management(outRas)
	outRas = None


in_path = out_path
out_path = landcover_path[:-4] + "_prp" + lc + ".tif"

#if not os.path.isfile(out_path) or not skip_existing:
#	#outRas = arcpy.sa.FocalStatistics(in_path,NbrRectangle(11,11,"CELL"),"MEAN","DATA")
#	outRas = arcpy.sa.FocalStatistics(in_path,NbrCircle(5,"CELL"),"MEAN","DATA")
#	#outRas.save(out_path)
#	outRas = arcpy.CopyRaster_management(outRas, out_path, "DEFAULTS", "0", "-999", "NONE", "NONE", "32_BIT_FLOAT")
#	outRas = arcpy.BuildPyramids_management(outRas)
#	outRas = None



##	Process distance-to rasters for landcovers of interest:
for lc in ['011', '040', '130', '140', '150', '160', '190', '200', '210', '230', '240', '250', 'BLT']:
	print("Distance:  " + lc)
	out_path = landcover_path[:-4] + "_dst" + lc + ".tif"

	if not os.path.isfile(out_path) or not skip_existing:
		lcRas = Raster(landcover_path[:-4] + "_cls" + lc + ".tif")
		tmpRas =  SetNull(lcRas != 1, 1)
		outRas = arcpy.sa.EucDistance(tmpRas)

		#outRas.save(out_path)
		outRas = arcpy.CopyRaster_management(outRas, out_path, "DEFAULTS", "0", "-999", "NONE", "NONE", "32_BIT_FLOAT")
		outRas = arcpy.BuildPyramids_management(outRas)
		outRas = None



if ("NPP" in dataset_folders):
	dataset_folders.remove("NPP")

	##	Clip and project our MODIS-derived NPP raster:
	dataset_folder = "NPP"
	print("PROCESSING:  " + dataset_folder)

	#dataset_name = "MOD17A3.A2010001.055.Npp_1km.tif"

	##	Instead of hardcoding it above we'll just pull in the only TIF or IMG
	##		file that should be in the directory:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.img" ) + glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.tif" ))
	if dataset_path:
		dataset_name = os.path.basename( dataset_path[0] )
	else:
		print("ERROR:  No NPP found yet the folder is present!  This indicates that there may be a problem with the NPP MODIS processing in the \"Data Preparation, R.r\" script or something else was overlooked!")
		exit()

	in_path = data_path + country + "/" + dataset_folder + "/" + dataset_name

	data_desc = arcpy.Describe(in_path)
	input_prj = data_desc.SpatialReference.exportToString()

	output_name = "prj.tif"

	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	if not os.path.isfile(out_path) or not skip_existing:
		##	NOTE: We use NEAREST as in the climate dataset because there are categorical values included in the dataset, and run it through a Con() to snap it to our land cover raster as the ProjectRaster_management() function does not honor it:
		arcpy.ProjectRaster_management(in_path,out_path,intermediate_prj,"NEAREST","100","#","#",input_prj)

		##	NOTE: So the ProjectRaster_management() function does seem to correctly
		##		apply the snapping and extent environment correctly, perhaps
		##		because we are not mosaicking at the same time?  In any case
		##		let's try running it without these, since with the nibbled BioClim
		##		data the Con() throws an error anyway:
		#outCon = Con(Raster(out_path), Raster(out_path), 0, "VALUE <> 0")

		##	Perform a couple more condtional statements to convert bare/built class
		##		values to 0 and Null for unspecified or error pixels:
		#outCon = Con(outCon, 0, outCon, "VALUE >= 65530 AND VALUE <= 65534")

		#outCon = Con(Raster(out_path), 0, Raster(out_path), "VALUE >= 65530 AND VALUE <= 65534")

		##	NOTE: The Con() statement fails more often than not, though
		##		inconsistently for some countries and for no reason that I
		##		can reliably troubleshoot.  Therefore we create a connection
		##		to the raster object and run the SetNull() functions first
		##		and this seems to more consistently finish.
		outCon = Raster(out_path)
		outCon = SetNull(outCon==65529, outCon)
		outCon = SetNull(outCon==65535, outCon)
		outCon = Con(outCon, 0, outCon, "VALUE >= 65530 AND VALUE <= 65534")

		##	TODO: When transitioning to Google Compute Engine we need to make sure that
		##		the NULL values get nibbled out with their nearest good value, like the
		##		WorldClim data...


		output_name = "npp.tif"
		out_path = tmp_path + "/" + output_name

		#outCon.save(out_path)
		outCon = arcpy.CopyRaster_management(outCon, out_path, "DEFAULTS", "0", "65535", "NONE", "NONE", "16_BIT_UNSIGNED")
		outCon = arcpy.BuildPyramids_management(outCon)
		outCon = None



if ("Lights" in dataset_folders):
	dataset_folders.remove("Lights")

	##	Clip and project our NASA Night Lights raster:
	dataset_folder = "Lights"
	print("PROCESSING:  " + dataset_folder)
	dataset_name = ""
	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder)

	##	If an image dataset exists in the folder we'll use it, otherwise we will
	##		pull the default dataset:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.img" ) + glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.tif" ))
	if dataset_path:
		dataset_name = os.path.basename( dataset_path[0] )
	else:
		dataset_name = ""

	if dataset_name == "":
		dataset_name = "NASA Night Lights.gdb/NASANightLightsMosaic"
		in_path = NightLights_path + dataset_name
	else:
		in_path = tmp_path + "/" + dataset_name

	##	Pull projection information for night lights dataset:
	data_desc = arcpy.Describe(in_path)
	input_prj = data_desc.SpatialReference.exportToString()


	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	mosaic_name = "mosaic.tif"
	mosaic_path = tmp_path + "/" + mosaic_name
	output_name = "mosaic_prj.tif"
	out_path = tmp_path + "/" + output_name


	if not os.path.isfile(out_path) or not skip_existing:
		##	 First we have to clear the environment spatial reference, extent and snapping:
		tmp_outputCoordinateSystem = arcpy.env.outputCoordinateSystem
		arcpy.env.outputCoordinateSystem = ""

		tmp_extent = arcpy.env.extent
		arcpy.env.extent = ""

		tmp_snapRaster = arcpy.env.snapRaster
		arcpy.env.snapRaster = ""


		##	Now that we can clip without any automatic reprojection occurring, we
		##		clip using our census buffer as our boundary extent:
		arcpy.Clip_management(in_path,"#",mosaic_path,buffer_path,"#","NONE")
		print("	FINISHED: Mosaic...")


		##	 Now we can restore the environment settings:
		##	 First we have to clear the environment spatial reference, extent and snapping:
		arcpy.env.outputCoordinateSystem = tmp_outputCoordinateSystem
		arcpy.env.extent = tmp_extent
		arcpy.env.snapRaster = tmp_snapRaster


		##	Next, we run the projection to project to our appropriate, output projection.  But note that the projection does't honor our output extent/snap environment settings, so we will run another conditional statement to fix this:

		##	NOTE: We use BILINEAR resampling here, not NEAREST as in the climate dataset, and run it through a Con() to snap it to our land cover raster as the ProjectRaster_management() function does not honor it:
		arcpy.ProjectRaster_management(mosaic_path, out_path, intermediate_prj,"BILINEAR","100","#","#",input_prj)

		##	NOTE: Because we don't specify using a particular band, the conditional statement only processes the first band (the red band of the three band composite) and so he output file only contains a single band representing the nonzero values of red.
		in_path = out_path

		output_name = "lights.tif"
		out_path = tmp_path + "/" + output_name

		arcpy.CopyRaster_management(in_path + "\\Band_1", out_path)
		print("	FINISHED: Clipping and projection...")



##	Now that we have rasters to snap to, we need to generate a zonal raster
##		from our administrative units based on the ADMINID, to be used for
##		calculation of zonal stats:
dataset_folder = "Census"
print("PROCESSING:  " + dataset_folder + " zonal data...")
dataset_name = "census_area.shp"

in_path = data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name

output_name = "census_zones.tif"
zonal_path = data_path + country + "/" + dataset_folder + "/Derived/" + output_name

if not os.path.isfile(zonal_path) or not skip_existing:
	##	NOTE: PolygonToRaster() does not honor the snap to raster environment
	##		settings, but FeatureToRaster_conversion() does, so we'll use it:
	arcpy.FeatureToRaster_conversion(in_path, "ADMINID", zonal_path, 100)


##	Convert our census blocks into points (located inside) for faster
##		zonal statistics processing in our R script:
in_path = in_path

output_name = "census_points.shp"
out_path = data_path + country + "/" + dataset_folder + "/Derived/" + output_name

if not os.path.isfile(out_path) or not skip_existing:
	arcpy.FeatureToPoint_management(in_path, out_path, "INSIDE")



if ("Roads" in dataset_folders):
	dataset_folders.remove("Roads")

	##	Extract, clip and project roads:
	dataset_folder = "Roads"
	print("PROCESSING:  " + dataset_folder)

	##	If a shapefile dataset exists in the folder we'll use it, otherwise we
	##		will pull the default dataset:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.shp" ))
	if dataset_path:
		dataset_name = os.path.basename( dataset_path[0] )
	else:
		dataset_name = ""

	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder)

	if dataset_name == "":
		dataset_name = "VMAP0.gdb/roadl"
		in_path = VMAP0_path + dataset_name
	else:
		in_path = tmp_path + "/" + dataset_name

	output_name = "roads.shp"
	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	if not os.path.isfile(out_path) or not skip_existing:
		##	Because ArcGIS 10.1 will process the extent before reprojection
		##		on the fly (which differs rom 9.3 and 10.0, which will reproject
		##		and use reprojected coordinates to compare to the provided extent)
		##		we must set the extent to None here, perform the clip (which
		##		reprojects on the fly) and then reset the default extent before
		##		proceeding:
		tmp_extent = arcpy.env.extent
		arcpy.env.extent = None
		arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
		arcpy.env.extent = tmp_extent

	process_point_area(dataset_folder)



if ("Rivers" in dataset_folders):
	dataset_folders.remove("Rivers")

	##	Extract, clip and project rivers/streams:
	dataset_folder = "Rivers"
	print("PROCESSING:  " + dataset_folder)

	##	If a shapefile dataset exists in the folder we'll use it, otherwise we
	##		will pull the default dataset:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.shp" ))
	if dataset_path:
		dataset_name = os.path.basename( dataset_path[0] )
	else:
		dataset_name = ""

	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder)

	if dataset_name == "":
		dataset_name = "VMAP0.gdb/watrcrsl"
		in_path = VMAP0_path + dataset_name
	else:
		in_path = data_path + country + "/" + dataset_folder + "/" + dataset_name

	output_name = "rivers.shp"
	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	if not os.path.isfile(out_path) or not skip_existing:
		##	See note above about clipping and process extent in ArcGIS 10.1:
		tmp_extent = arcpy.env.extent
		arcpy.env.extent = None
		arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
		arcpy.env.extent = tmp_extent

	process_point_area(dataset_folder)



if ("Waterbodies" in dataset_folders):
	dataset_folders.remove("Waterbodies")

	##	Extract, clip and project water bodies:
	dataset_folder = "Waterbodies"
	print("PROCESSING:  " + dataset_folder)

	##	If a shapefile dataset exists in the folder we'll use it, otherwise we
	##		will pull the default dataset:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.shp" ))
	if dataset_path:
		dataset_name = os.path.basename( dataset_path[0] )
	else:
		dataset_name = ""

	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder)

	if dataset_name == "":
		dataset_name = "VMAP0.gdb/inwatera"
		in_path = VMAP0_path + dataset_name
	else:
		in_path = data_path + country + "/" + dataset_folder + "/" + dataset_name

	output_name = "waterbodies.shp"
	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	if not os.path.isfile(out_path) or not skip_existing:
		##	See note above about clipping and process extent in ArcGIS 10.1:
		tmp_extent = arcpy.env.extent
		arcpy.env.extent = None
		arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
		arcpy.env.extent = tmp_extent


	process_point_area(dataset_folder)



if ("Populated" in dataset_folders):
	dataset_folders.remove("Populated")

	##	Process populated areas dataset(s):
	dataset_folder = "Populated"
	print("PROCESSING:  " + dataset_folder)

	##	If a shapefile dataset exists in the folder we'll use it, otherwise we
	##		will pull the default dataset:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.shp" ))
	if dataset_path:
		dataset_name = os.path.basename( dataset_path[0] )
	else:
		dataset_name = ""

	if dataset_name == "":
		tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")

		output_name = "populated.shp"
		out_path = tmp_path + "/" + output_name

		if not os.path.isfile(out_path) or not skip_existing:

			##	Set an empty string to hold the list of area-based datasets to merge at the end of this:
			merge_paths = ""


			##	Extract, clip and project miscellaneous populated point locations and built areas:

			dataset_name = "VMAP0.gdb/builtupa"
			in_path = VMAP0_path + dataset_name
			output_name = "populated_builtupa.shp"
			out_path = tmp_path + "/" + output_name

			##	See note above about clipping and process extent in ArcGIS 10.1:
			tmp_extent = arcpy.env.extent
			arcpy.env.extent = None
			arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
			arcpy.env.extent = tmp_extent

			merge_paths += "'" + out_path + "';"


			##	Removed because it doesn't exist in the South America/Africa VMAP0
			##		data layers:
			#dataset_name = "VMAP0.gdb/mispopa"
			#in_path = VMAP0_path + dataset_name
			#output_name = "populated_mispopa.shp"
			#out_path = tmp_path + "/" + output_name
			#
			###	See note above about clipping and process extent in ArcGIS 10.1:
			#tmp_extent = arcpy.env.extent
			#arcpy.env.extent = None
			#arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
			#arcpy.env.extent = tmp_extent
			#
			#merge_paths += "'" + out_path + "';"


			##	Removed because it doesn't exist in the South America/Africa VMAP0
			##		data layers:
			#dataset_name = "VMAP0.gdb/builtupp"
			#in_path = VMAP0_path + dataset_name
			#output_name = "populated_builtupp.shp"
			#out_path = tmp_path + "/" + output_name
			#
			###	See note above about clipping and process extent in ArcGIS 10.1:
			#tmp_extent = arcpy.env.extent
			#arcpy.env.extent = None
			#arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
			#arcpy.env.extent = tmp_extent


			dataset_name = "VMAP0.gdb/mispopp"
			in_path = VMAP0_path + dataset_name
			output_name = "populated_mispopp.shp"
			out_path = tmp_path + "/" + output_name

			##	See note above about clipping and process extent in ArcGIS 10.1:
			tmp_extent = arcpy.env.extent
			arcpy.env.extent = None
			arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
			arcpy.env.extent = tmp_extent


			##	Buffer built and miscellaneous populated point locations and add
			##		to the merge list:

			##	Removed because it doesn't exist in the South America/Africa VMAP0
			##		data layers:
			#dataset_name = "populated_builtupp.shp"
			#in_path = tmp_path + "/" + dataset_name
			#output_name = "populated_builtupp_buff.shp"
			#out_path = tmp_path + "/" + output_name
			#arcpy.Buffer_analysis(in_path, out_path,"100 Meters","FULL","ROUND","NONE","#")
			#merge_paths += "'" + out_path + "';"

			dataset_name = "populated_mispopp.shp"
			in_path = tmp_path + "/" + dataset_name
			output_name = "populated_mispopp_buff.shp"
			out_path = tmp_path + "/" + output_name
			arcpy.Buffer_analysis(in_path, out_path,"100 Meters","FULL","ROUND","NONE","#")
			merge_paths += "'" + out_path + "';"


			##	Merge populated areas:
			output_name = "populated.shp"
			out_path = tmp_path + "/" + output_name
			arcpy.Merge_management(merge_paths, out_path,"#")


	else:
		tmp_path = ensure_dir(data_path + country + "/" + dataset_folder)
		in_path = tmp_path + "/" + dataset_name

		output_name = "populated.shp"
		tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
		out_path = tmp_path + "/" + output_name

		if not os.path.isfile(out_path) or not skip_existing:
			##	See note above about clipping and process extent in ArcGIS 10.1:
			tmp_extent = arcpy.env.extent
			arcpy.env.extent = None
			arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
			arcpy.env.extent = tmp_extent


	process_point_area(dataset_folder)



if ("Protected" in dataset_folders):
	dataset_folders.remove("Protected")

	##	Extract by clip and project protected areas:
	dataset_folder = "Protected"
	print("PROCESSING:  " + dataset_folder)

	##	If a shapefile dataset exists in the folder we'll use it, otherwise we
	##		will pull the default dataset:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.shp" ))
	if dataset_path:
		dataset_name = os.path.basename( dataset_path[0] )
	else:
		dataset_name = ""

	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder)

	if dataset_name == "":
		dataset_name = "WDPAfgdb_Sept2012.gdb/WDPApoly_September2012"
		in_path = WDPA_path + dataset_name
	else:
		in_path = tmp_path + "/" + dataset_name

	output_name = "protected.shp"
	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	if not os.path.isfile(out_path) or not skip_existing:
		##	See note above about clipping and process extent in ArcGIS 10.1:
		tmp_extent = arcpy.env.extent
		arcpy.env.extent = None
		arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
		arcpy.env.extent = tmp_extent

	process_point_area(dataset_folder)



if ("Urban" in dataset_folders):
	dataset_folders.remove("Urban")

	##	Extract by clip and project urban extent areas (Schneider):
	dataset_folder = "Urban"
	print("PROCESSING:  " + dataset_folder)


	##	Check to see if instead of a shapefile urban extent we have a raster
	##		for Urban data:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.img" ) + glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.tif" ))

	if dataset_path:
		##	If we have a raster, we're just going to add this back in to
		##		process as a non-default data set and skip the rest of the
		##		processing here:
		dataset_folders.append("Urban")

	else:
		##	If a shapefile dataset exists in the folder we'll use it,
		##		otherwise we will pull the default dataset:
		dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.shp" ))
		if dataset_path:
			dataset_name = os.path.basename( dataset_path[0] )
		else:
				dataset_name = ""

		tmp_path = ensure_dir(data_path + country + "/" + dataset_folder)

		if dataset_name == "":
			dataset_name = (glob.glob( UrbExtents_path + "/" + "*.shp" ))
			in_path = dataset_name[0]
		else:
			in_path = tmp_path + "/" + dataset_name

		output_name = "urban.shp"
		tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
		out_path = tmp_path + "/" + output_name

		if not os.path.isfile(out_path) or not skip_existing:
			##	See note above about clipping and process extent in ArcGIS 10.1:
			tmp_extent = arcpy.env.extent
			arcpy.env.extent = None
			arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
			arcpy.env.extent = tmp_extent

		process_point_area(dataset_folder)



if ("Elevation" in dataset_folders):
	dataset_folders.remove("Elevation")

	##	Reproject and clip our raster catalog/mosaic of HydroSHEDS tiles:
	dataset_folder = "Elevation"
	print("PROCESSING:  " + dataset_folder)

	##	Instead of hardcoding it above we'll look for a custom TIF or IMG
	##		file that should be in the directory:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.img" ) + glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.tif" ))
	if dataset_path:
		dataset_name = os.path.basename( dataset_path[0] )
	else:
		dataset_name = ""

	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder)

	if dataset_name == "":
		dataset_name = "Void-Filled DEM.gdb/VoidFilledDEMMosaic"
		in_path = HydroSHEDS_path + dataset_name
	else:
		in_path = tmp_path + "/" + dataset_name

	##	Pull projection information for base elevation dataset:
	data_desc = arcpy.Describe(in_path)
	input_prj = data_desc.SpatialReference.exportToString()


	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	mosaic_name = "mosaic.tif"
	mosaic_path = tmp_path + "/" + mosaic_name
	output_name = "elevation.tif"
	out_path = tmp_path + "/" + output_name


	if not os.path.isfile(out_path) or not skip_existing:
		##	 First we have to clear the environment spatial reference, extent and snapping:
		tmp_outputCoordinateSystem = arcpy.env.outputCoordinateSystem
		arcpy.env.outputCoordinateSystem = ""

		tmp_extent = arcpy.env.extent
		arcpy.env.extent = ""

		tmp_snapRaster = arcpy.env.snapRaster
		arcpy.env.snapRaster = ""


		##	Now that we can clip without any automatic reprojection occurring, we
		##		clip using our census buffer as our boundary extent:
		arcpy.Clip_management(in_path,"#",mosaic_path,buffer_path,"#","NONE")
		print("	FINISHED: Mosaic...")


		##	 Now we can restore the environment settings:
		##	 First we have to clear the environment spatial reference, extent and snapping:
		arcpy.env.outputCoordinateSystem = tmp_outputCoordinateSystem
		arcpy.env.extent = tmp_extent
		arcpy.env.snapRaster = tmp_snapRaster


		##	Next, we run the projection to project to our appropriate, output projection.  But note that the projection does't honor our output extent/snap environment settings, so we will run another conditional statement to fix this:

		##	NOTE: We use BILINEAR resampling here, not NEAREST as in the climate dataset, and run it through a Con() to snap it to our land cover raster as the ProjectRaster_management() function does not honor it:
		arcpy.ProjectRaster_management(mosaic_path, out_path, intermediate_prj,"BILINEAR","100","#","#",input_prj)
		##	NOTE: See above about the Con() statement and its removal...
		#outCon = Con(Raster(out_path), Raster(out_path), 0, "VALUE <> 0")
		#outCon.save(out_path)
		#outCon = None
		print("	FINISHED: Projection...")


		##	Now create our slope data from this file:
		outRas = Slope(out_path)

		output_name = "elevation_slope.tif"
		out_path = tmp_path + "/" + output_name

		#outRas.save(out_path)
		outRas = arcpy.CopyRaster_management(outRas, out_path, "DEFAULTS", "0", "-999", "NONE", "NONE", "32_BIT_FLOAT")
		outRas = arcpy.BuildPyramids_management(outRas)
		outRas = None
		print("	FINISHED: Slope...")



if ("Temp" in dataset_folders):
	dataset_folders.remove("Temp")

	##	Reproject and clip our WorldClim temperature data to match our landcover:
	dataset_folder = "Temp"
	print("PROCESSING:  " + dataset_folder)

	#	Instead of hardcoding it above we'll look for a custom TIF or IMG
	#		file that should be in the directory:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.img" ) + glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.tif" ))
	if dataset_path:
		dataset_name = os.path.basename( dataset_path[0] )
	else:
		dataset_name = ""

	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder)

	if dataset_name == "":
		dataset_name = "bio_1_nibble.img"
		in_path = WorldClim_path + dataset_name
	else:
		in_path = data_path + country + "/" + dataset_folder + "/" + dataset_name

	##	Pull projection information for base dataset:
	data_desc = arcpy.Describe(in_path)
	input_prj = data_desc.SpatialReference.exportToString()

	output_name = "temp.tif"
	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	if not os.path.isfile(out_path) or not skip_existing:
		##	 First we have to clear the environment spatial reference, extent and snapping:
		tmp_outputCoordinateSystem = arcpy.env.outputCoordinateSystem
		arcpy.env.outputCoordinateSystem = ""

		tmp_extent = arcpy.env.extent
		arcpy.env.extent = ""

		tmp_snapRaster = arcpy.env.snapRaster
		arcpy.env.snapRaster = ""


		##	Now that we can clip without any automatic reprojection occurring, we
		##		clip using our census buffer as our boundary extent:
		mosaic_name = "mosaic.tif"
		mosaic_path = tmp_path + "/" + mosaic_name

		arcpy.Clip_management(in_path,"#",mosaic_path,buffer_path,"#","NONE")
		print("	FINISHED: Mosaic...")


		##	 Now we can restore the environment settings:
		##	 First we have to clear the environment spatial reference, extent and snapping:
		arcpy.env.outputCoordinateSystem = tmp_outputCoordinateSystem
		arcpy.env.extent = tmp_extent
		arcpy.env.snapRaster = tmp_snapRaster


		##	Next, we run the projection to project to our appropriate, output projection.  But note that the projection does't honor our output extent/snap environment settings, so we will run another conditional statement to fix this:

		##	NOTE: We use BILINEAR as in the elevation dataset, and need to use the Con() function in order to snap it to our land cover, because ProjectRaster_management() does not honor it:
		arcpy.ProjectRaster_management(mosaic_path, out_path, intermediate_prj,"BILINEAR","100","#","#",input_prj)
		##	NOTE: These Con() statements were buggy with the nibbled BioClim data to
		##		begin with.  They wouldn't even run in ArcMap, causing Python
		##		to crash immediately upon running.
		#outCon = Con(Raster(out_path), Raster(out_path), 0, "VALUE <> 0")
		#outCon.save(out_path)
		#outCon = None



if ("Precip" in dataset_folders):
	dataset_folders.remove("Precip")

	##	Reproject and clip our WorldClim precipitation data to match our landcover:
	dataset_folder = "Precip"
	print("PROCESSING:  " + dataset_folder)

	#	Instead of hardcoding it above we'll look for a custom TIF or IMG
	#		file that should be in the directory:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.img" ) + glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.tif" ))
	if dataset_path:
		dataset_name = os.path.basename( dataset_path[0] )
	else:
		dataset_name = ""

	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder)

	if dataset_name == "":
		dataset_name = "bio_12_nibble.img"
		in_path = WorldClim_path + dataset_name
	else:
		in_path = data_path + country + "/" + dataset_folder + "/" + dataset_name

	##	Pull projection information for base dataset:
	data_desc = arcpy.Describe(in_path)
	input_prj = data_desc.SpatialReference.exportToString()

	output_name = "precip.tif"
	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	if not os.path.isfile(out_path) or not skip_existing:
		##	 First we have to clear the environment spatial reference, extent and snapping:
		tmp_outputCoordinateSystem = arcpy.env.outputCoordinateSystem
		arcpy.env.outputCoordinateSystem = ""

		tmp_extent = arcpy.env.extent
		arcpy.env.extent = ""

		tmp_snapRaster = arcpy.env.snapRaster
		arcpy.env.snapRaster = ""


		##	Now that we can clip without any automatic reprojection occurring, we
		##		clip using our census buffer as our boundary extent:
		mosaic_name = "mosaic.tif"
		mosaic_path = tmp_path + "/" + mosaic_name

		arcpy.Clip_management(in_path,"#",mosaic_path,buffer_path,"#","NONE")
		print("	FINISHED: Mosaic...")


		##	 Now we can restore the environment settings:
		##	 First we have to clear the environment spatial reference, extent and snapping:
		arcpy.env.outputCoordinateSystem = tmp_outputCoordinateSystem
		arcpy.env.extent = tmp_extent
		arcpy.env.snapRaster = tmp_snapRaster


		##	Next, we run the projection to project to our appropriate, output projection.  But note that the projection does't honor our output extent/snap environment settings, so we will run another conditional statement to fix this:

		##	NOTE: We use BILINEAR as in the elevation dataset, and need to use the Con() function in order to snap it to our land cover, because ProjectRaster_management() does not honor it:
		arcpy.ProjectRaster_management(mosaic_path, out_path, intermediate_prj,"BILINEAR","100","#","#",input_prj)
		##	NOTE: These Con() statements were buggy with the nibbled BioClim data to
		##		begin with.  They wouldn't even run in ArcMap, causing Python
		##		to crash immediately upon running.
		#outCon = Con(Raster(out_path), Raster(out_path), 0, "VALUE <> 0")
		#outCon.save(out_path)
		#outCon = None


##	END: Data pre-processing for default datasets
#####



#####
##	BEGIN: Data pre-processing for non-default datasets

##	Any additional ancillary datasets will remain in the dataset_folders
##		array, and we are going to process them according to their detected
##		data type:

for dataset_folder in dataset_folders:

	##	Check to see if this is a folder to skip:
	if not dataset_folder.startswith("!"):

		dataset_name = dataset_folder.lower()

		print("PROCESSING:  " + dataset_folder + " Non-Default Data")


		##	Determine the type of data present (raster or shapefile):


		##	If a shapefile dataset exists in the folder we'll use it, otherwise we
		##		will exit with an error:
		dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.shp" ))
		if dataset_path:
			dataset_orig_name = os.path.basename( dataset_path[0] )
			dataset_type = "shapefile"

		else:
			dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.img" ) + glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.tif" ))
			dataset_type = "raster"

			if dataset_path:
				dataset_orig_name = os.path.basename( dataset_path[0] )

			else:
				print("ERROR:  No data file (*.img, *.tif, or *.shp) found in " + dataset_folder + "!")
				exit()


		tmp_path = ensure_dir(data_path + country + "/" + dataset_folder)
		in_path = tmp_path + "/" + dataset_orig_name


		if dataset_type == "shapefile":
			output_name = dataset_name + ".shp"
			tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
			out_path = tmp_path + "/" + output_name

			if not os.path.isfile(out_path) or not skip_existing:
				##	See note above about clipping and process extent in ArcGIS 10.1:
				tmp_extent = arcpy.env.extent
				arcpy.env.extent = None
				arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
				arcpy.env.extent = tmp_extent

			process_point_area(dataset_folder)

		else:
			##	Process as a raster:
			dataset_path = data_path + country + "/" + dataset_folder + "/" + dataset_orig_name

			##	Determine whether the raster is binary or continuous:

			##	First see if this is a floating point raster:
			value_type_result = arcpy.GetRasterProperties_management(dataset_path, "VALUETYPE")
			value_type = value_type_result.getOutput(0)

			if value_type == "9" or value_type == "10":
				##	If it is floating point or double then process as continuous:
				process_raster_continuous(dataset_folder)
			else:
				##	If not, then check to see if it has more than 2 unique values:
				unique_value_count_result = arcpy.GetRasterProperties_management(dataset_path, "UNIQUEVALUECOUNT")
				unique_value_count = unique_value_count_result.getOutput(0)

				if int(unique_value_count) <= 3:
					process_raster_binary(dataset_folder)
				else:
					process_raster_continuous(dataset_folder)


##	END: Data pre-processing for non-default datasets
#####



#####
##	Preprocessing Notes for randomForest model:
##		1)	Land Cover and Other Areal Information (e.g. Protected Area)
##			-	Broken down into proportion of each land cover class by square
##				kilometer using a moving window of 1000 m x 1000 m for each 100 m
##				pixel
##			-	Typical data:
##					-	Land cover class:
##						11 	= Cultivated Terrestrial Areas and Managed Lands
##						40  = Natural and Semi-natural Terrestrial Vegetation ? Woody / Trees
##						130 = Natural and Seminatural Terrestrial Vegetation ? Shrubs
##						140 = Natural and Seminatural Terrestrial Vegetation ? Herbaceous
##						150 = Natural and Semi-natural Terrestrial Vegetation
##						160 = Natural and Seminatural Aquatic Vegetation
##						190 = Urban area
##						200 = Bare areas
##						210 = Water bodies
##						230 = No data, cloud/shadow
##						240	= Rural settlement
##						250	= Industrial area
##						BLT = Combined 190 and 240 built areas
##					-	Protected Area
##						-	Default: World Database on Protected Areas, UN
##					-	Urban Extents
##						-	Default: Classified from MODIS-derived data (Schneider, et al., UN)
##					-	Waterbodies
##						-	Default: vmaplv0/sasaus:hydro/inwatera
##					-	Populated Areas, Merged (points are buffered by 100 m)
##						-	Default:
##								vmaplv0/sasaus:pop/builtupp (buffered to 100 m)
##								vmaplv0/sasaus:pop/builtupa
##								vmaplv0/sasaus:pop/mispopp (buffered to 100 m)
##								(vmaplv0/sasaus:pop/mispopa) [REMOVED]
##									Because	it's not present in the South	America and
##										Africa dataset
##						-	But may include OSM, health care centers, and any variety
##							of ancillary data collected for each country...
##			-	These pixels are averaged per census block for parameterizing
##				the random forest model
##
##		2)	Distance-to Land Cover Class, Area Polygon, or Other Feature
##			-	A distance-to raster is created based on a binary mask of several
##				different types of raster-, polygon- or line- based data
##			-	Typical data:
##				-	Land cover classes:
##						11 	= Cultivated Terrestrial Areas and Managed Lands
##						40  = Natural and Semi-natural Terrestrial Vegetation ? Woody / Trees
##						130 = Natural and Seminatural Terrestrial Vegetation ? Shrubs
##						140 = Natural and Seminatural Terrestrial Vegetation ? Herbaceous
##						150 = Natural and Semi-natural Terrestrial Vegetation
##						160 = Natural and Seminatural Aquatic Vegetation
##						190 = Urban area
##						200 = Bare areas
##						210 = Water bodies
##						230 = No data, cloud/shadow
##						240	= Rural settlement
##						250	= Industrial area
##						BLT = Combined 190 and 240 built areas
##				-	Rivers
##					-	Default: vmaplv0/sasaus:hydro/watrcrsl
##				-	Roads
##					-	Default: vmaplv0/sasaus:trans/roadl
##				-	Waterbodies (see above)
##				-	Populated Areas, Merged (see above)
##				-	Protected Areas (see above)
##				-	Urban Extent Areas (see above)
##
##		3)	Mean value per administrative unit:
##			- Typical data:
##				-	Elevation
##				-	Slope
##				-	V1 Afri/AsiaPop Prediction Density Weightings
##					-	Default: gridx
##				-	MODIS-derived NPP (MOD17A3)
##				-	NASA Night Lights, derived from Suomi NPP VIIRS sensed mosaics
#####
