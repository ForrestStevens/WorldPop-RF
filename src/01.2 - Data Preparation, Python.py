#####
##	BEGIN:	Set per-country configuration options

country = "NGA"
#country = "KHM_002"
#country = "KHM_002_Anc"


##	Configure project and default data folders:
root_path = "D:/Documents/Graduate School/Research/Population/Data/"
#root_path = "C:/Users/forrest/Research/Population/Data/"
#root_path = "D:/Research/Population/Data/"


##	NOTE: This needs to be adjusted depending on whether we are in
##		the Asia/Australia region or whether we are in the Africa/America
##		region.  If you are processing a large country that does not fall
##		*completely* within one of he VMAP0 data regions then you should
##		either merge the datasets by hand or provide alternatives to the 
##		following default datasets:  "Roads", "Rivers", "Populated"
#NGA_path = root_path + "GIS/NGA/VMAP0/v0sas_5/vmaplv0/sasaus:"
NGA_path = root_path + "GIS/NGA/VMAP0/v0soa_5/vmaplv0/soamafr:"
#NGA_path = root_path + "GIS/NGA/VMAP0/v0noa_5/vmaplv0/noamer:"
#NGA_path = root_path + "GIS/NGA/VMAP0/v0eur_5/vmaplv0/eurnasia:"


##	Output projection is configured for each country based on the most
##		appropriate output for distance/area considerations:
##	For KHM, VNM:
#intermediate_prj = "PROJCS['WGS_1984_UTM_Zone_48N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',105.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]"
##	For KEN:
#intermediate_prj = "PROJCS['WGS_1984_UTM_Zone_37N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',39.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]"
##      For NGA:
intermediate_prj = "PROJCS['WGS_1984_UTM_Zone_32N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',9.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0],AUTHORITY['EPSG',32632]]"


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

data_path = root_path + "RF/data/"
output_path = root_path + "RF/output/"

#Peel_path = root_path + "GIS/Peel Climate Classification/Recoded/"
WorldClim_path = root_path + "WorldClim/BioClim/"

HydroSHEDS_path = root_path + "HydroSHEDS/Void-Filled DEM/"
WDPA_path = root_path + "GIS/Protected Areas/"

#NightLights_path = root_path + "NASA Night Lights/"
NightLights_path = root_path + "NASA Night Lights, Unprocessed/"

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

##	Should we overwrite any already existing derived data:
overwrite = True
arcpy.env.overwriteOutput = overwrite


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
	outPath = data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + "_dst" +".tif"
	
	if not os.path.isfile(outPath) or not skip_existing:
		outRas = arcpy.sa.EucDistance(data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + ".shp")
		outRas.save(outPath)


def process_point_area(dataset_folder):
	##	Process select feature data into binary masks, distance-to
	##		and proportions per	81 hectares (11 x 11 cell moving window, 
	##		specified as a NbrCircle with radius 5):
	
	dataset_folder = dataset_folder
	dataset_name = dataset_folder.lower()

	print("PROCESSING:  " + dataset_folder + " Features")

	##	Use FID to create output raster from our features:
	tmpRas = data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + "_FID.tif"

	if not os.path.isfile(tmpRas) or not skip_existing:
		##	See note below about PolygonToRaster and it not respecting the 
		##		snap-to raster environment... Same here with all FeatureToRaster
		##		calls...
		arcpy.FeatureToRaster_conversion(data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + ".shp", "FID", tmpRas, 100)

		##	Create distance-to raster:
		outRas = arcpy.sa.EucDistance(tmpRas)
		outRas.save(data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + "_dst.tif")

		##	Save the binary classification raster:
		outRas =  IsNull(Raster(tmpRas)) == 0
		outRas.save(data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + "_cls.tif")

		##	Calculate proportion of cover:
		#outRas = arcpy.sa.FocalStatistics(outRas,NbrRectangle(11,11,"CELL"),"MEAN","DATA")
		outRas = arcpy.sa.FocalStatistics(outRas,NbrCircle(5,"CELL"),"MEAN","DATA")
		outRas.save(data_path + country + "/" + dataset_folder + "/Derived/" + dataset_name + "_prp.tif")


def process_raster_binary(dataset_folder):
	dataset_folder = dataset_folder
	dataset_name = dataset_folder.lower()

	#	Clip and project our built raster:
	print("PROCESSING:  " + dataset_folder)

	#	Instead of hardcoding it above we'll look for a custom TIF or IMG 
	#		file that should be in the directory:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.img" ) + glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.tif" ))
	if dataset_path:
		dataset_name = os.path.basename( dataset_path[0] )
	else:
		dataset_name = ""

	in_path = data_path + country + "/" + dataset_folder + "/" + dataset_name

	output_name = dataset_name + "_cls.tif"

	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	data_desc = arcpy.Describe(in_path)
	input_prj = data_desc.SpatialReference.exportToString()

	
	if not os.path.isfile(out_path) or not skip_existing:
		arcpy.ProjectRaster_management(in_path,out_path,intermediate_prj,"NEAREST","100","#","#",input_prj)

		print("Proportion:  " + dataset_folder)
		#outRas = arcpy.sa.FocalStatistics(out_path,NbrRectangle(11,11,"CELL"),"MEAN","DATA")
		outRas = arcpy.sa.FocalStatistics(out_path,NbrCircle(5,"CELL"),"MEAN","DATA")

		out_path = out_path[:-8] + "_prp.tif"
		outRas.save(out_path)


		print("Distance:  " + dataset_folder)
		lcRas = Raster(out_path[:-8] + "_cls.tif")
		tmpRas =  SetNull(lcRas != 1, 1)
		outRas = arcpy.sa.EucDistance(tmpRas)

		out_path = out_path[:-8] + "_dst.tif"
		outRas.save(out_path)


def process_raster_continuous(dataset_folder):
	dataset_folder = dataset_folder
	dataset_name = dataset_folder.lower()

	#	Clip and project our built raster:
	print("PROCESSING:  " + dataset_folder)

	#	Instead of hardcoding it above we'll look for a custom TIF or IMG 
	#		file that should be in the directory:
	dataset_path = (glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.img" ) + glob.glob( data_path + country + "/" + dataset_folder + "/" + "*.tif" ))
	if dataset_path:
		dataset_name = os.path.basename( dataset_path[0] )
	else:
		dataset_name = ""

	in_path = data_path + country + "/" + dataset_folder + "/" + dataset_name

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


##	Now that we have a tmp_path specified, load the JSON objects that
##		specify the metadata and covariates specific to each dataset we
##		need to process:
metadata_json = open(tmp_path + "metadata.json")
covariates_json = open(tmp_path + "covariates.json")

metadata = json.load(metadata_json)
covariates = json.load(covariates_json)


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
		outRas.save(out_path)

	in_path = out_path
	out_path = landcover_path[:-4] + "_prp" + lc + ".tif"
	if not os.path.isfile(out_path) or not skip_existing:
		#outRas = arcpy.sa.FocalStatistics(in_path,NbrRectangle(11,11,"CELL"),"MEAN","DATA")
		outRas = arcpy.sa.FocalStatistics(in_path,NbrCircle(5,"CELL"),"MEAN","DATA")
		outRas.save(out_path)


##	Process land cover combinations of interest:

##	Combined built class (land cover 190 and 240):
lc = "BLT"
print(lc)
out_path = landcover_path[:-4] + "_cls" + lc + ".tif"
if not os.path.isfile(out_path) or not skip_existing:
	outRas =  (Raster(landcover_path) == 190) + (Raster(landcover_path) == 240)
	outRas.save(out_path)


in_path = out_path
out_path = landcover_path[:-4] + "_prp" + lc + ".tif"

if not os.path.isfile(out_path) or not skip_existing:
	#outRas = arcpy.sa.FocalStatistics(in_path,NbrRectangle(11,11,"CELL"),"MEAN","DATA")
	outRas = arcpy.sa.FocalStatistics(in_path,NbrCircle(5,"CELL"),"MEAN","DATA")
	outRas.save(out_path)



##	Process distance-to rasters for landcovers of interest:
for lc in ['011', '040', '130', '140', '150', '160', '190', '200', '210', '230', '240', '250', 'BLT']:
	print("Distance:  " + lc)
	out_path = landcover_path[:-4] + "_dst" + lc + ".tif"
	
	if not os.path.isfile(out_path) or not skip_existing:
		lcRas = Raster(landcover_path[:-4] + "_cls" + lc + ".tif")
		tmpRas =  SetNull(lcRas != 1, 1)
		outRas = arcpy.sa.EucDistance(tmpRas)

		outRas.save(out_path)



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

		##      NOTE: The Con() statement fails more often than not, though
		##              inconsistently for some countries and for no reason that I
		##              can reliably troubleshoot.  Therefore we create a connection
		##              to the raster object and run the SetNull() functions first
		##              and this seems to more consistently finish.
		outCon = Raster(out_path)
		outCon = SetNull(outCon==65529, outCon)
		outCon = SetNull(outCon==65535, outCon)
		outCon = Con(outCon, 0, outCon, "VALUE >= 65530 AND VALUE <= 65534")


		output_name = "npp.tif"
		out_path = tmp_path + "/" + output_name

		outCon.save(out_path)



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
		dataset_name = "trans/roadl"
		in_path = NGA_path + dataset_name
	else:
		in_path = tmp_path + "/" + dataset_name

	output_name = "roads.shp"
	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	if not os.path.isfile(out_path) or not skip_existing:
		arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")

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
		dataset_name = "hydro/watrcrsl"
		in_path = NGA_path + dataset_name
	else:
		in_path = data_path + country + "/" + dataset_folder + "/" + dataset_name

	output_name = "rivers.shp"
	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	if not os.path.isfile(out_path) or not skip_existing:
		arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
	
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
		dataset_name = "hydro/inwatera"
		in_path = NGA_path + dataset_name
	else:
		in_path = data_path + country + "/" + dataset_folder + "/" + dataset_name

	output_name = "waterbodies.shp"
	tmp_path = ensure_dir(data_path + country + "/" + dataset_folder + "/Derived")
	out_path = tmp_path + "/" + output_name

	if not os.path.isfile(out_path) or not skip_existing:
		arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
	
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
			
			dataset_name = "pop/builtupa"
			in_path = NGA_path + dataset_name
			output_name = "populated_builtupa.shp"
			out_path = tmp_path + "/" + output_name
			arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
			merge_paths += "'" + out_path + "';"

			##	Removed because it doesn't exist in the South America/Africa VMAP0
			##		data layers:
			#dataset_name = "pop/mispopa"
			#in_path = NGA_path + dataset_name
			#output_name = "populated_mispopa.shp"
			#out_path = tmp_path + "/" + output_name
			#arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
			#merge_paths += "'" + out_path + "';"

			
			##	Removed because it doesn't exist in the South America/Africa VMAP0
			##		data layers:
			#dataset_name = "pop/builtupp"
			#in_path = NGA_path + dataset_name
			#output_name = "populated_builtupp.shp"
			#out_path = tmp_path + "/" + output_name
			#arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")

			dataset_name = "pop/mispopp"
			in_path = NGA_path + dataset_name
			output_name = "populated_mispopp.shp"
			out_path = tmp_path + "/" + output_name
			arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")

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
			arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")
	
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
		arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")

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
		print("	FINISHED: Projection...")


		##	Now create our slope data from this file:
		outRas = Slope(out_path)

		output_name = "elevation_slope.tif"
		out_path = tmp_path + "/" + output_name

		outRas.save(out_path)
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
				arcpy.Clip_analysis(in_path, buffer_path, out_path,"#")

			process_point_area(dataset_folder)
		
		else:
			##	Process as a raster:

			##	Determine whether the raster is binary or continuous:
			unique_value_count_result = arcpy.GetRasterProperties_management(dataset_orig_name, "UNIQUEVALUECOUNT")
			unique_value_count = unique_value_count_result.getOutput(0) 

			if unique_value_count <= 3:
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
##						40  = Natural and Semi-natural Terrestrial Vegetation – Woody / Trees
##						130 = Natural and Seminatural Terrestrial Vegetation – Shrubs
##						140 = Natural and Seminatural Terrestrial Vegetation – Herbaceous
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
##						40  = Natural and Semi-natural Terrestrial Vegetation – Woody / Trees
##						130 = Natural and Seminatural Terrestrial Vegetation – Shrubs
##						140 = Natural and Seminatural Terrestrial Vegetation – Herbaceous
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
