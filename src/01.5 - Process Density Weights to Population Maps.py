#####
##	NOTE: In order for this script to run as both a stand-alone Python
##		geoprocessing script, and within ArcCatalog/ArcMap via "Load..."
##		in the Python interpreter window I added arcpy.MakeFeatureLayer...
##		calls after the creation of temporary datasets.	This unfortunately
##		breaks when you're running in ArcMap which by default adds the
##		temporary datasets to the current document.	To fix this I added
##		checks at each step to see if that exists in order for it to process.
#####

#####
## Purpose of script: Project final population map (popmap) to 2010 and 
##		2015 estimates based on the 2011 Revision of World Urbanization 
##		Prospects (http://esa.un.org/unpd/wup/index.htm) and also adjusted 
##		to U.N. total population numbers for each country.

##	Population/census data shapefile should be updated with columns
##		'ADMINID', 'ADMINPOP' and 'ISO'
##	(NOTE: takes place of admin.shp in Benin example)

##	GRUMP shapefile is set below, e.g. af_as_lac_urban.shp, a subset of
##		GRUMP world shapefile with columns 'ISOURBID' and 'POP' used in code
##		(already in shapefile)

##	Urban and rural growth rates to adjust to 2010 and also 2015
##		(growth rates excel file)

##	Other notes:
##		Check in initial population census shapefile for a RECNO and POP fields - if there are, delete
##	Geometry of census shp should be checked. if issues, topology may need to be corrected/adjusted (ex. VNM)
#####



#####
##	BEGIN:	Set per-country configuration options

country = "NGA"
#country = "KHM_002"
#country = "KHM_002_Anc"


##	Configure project and default data folders:
root_path = "D:/Documents/Graduate School/Research/Population/Data/"
#root_path = "C:/Users/forrest/Research/Population/Data/"
#root_path = "D:/Research/Population/Data/"

data_path = root_path + "RF/data/" + country + "/"
output_path = root_path + "RF/output/" + country + "/"
tmp_path = output_path + "tmp/"


##	TODO: Eventually this can be set in the Metadata.r file and pulled
##		via JSON but for now we will set it here, you must make
##              sure it matches the versioning in the Metadata.r file:

##	Versioning information:
version = "2b"


##	Country specific population-specific variables:
##		Growth rates are estimated for 2010 and 2015. In this example, they
##		are listed as 1 because country population census is 2010 and estimated
##		for 2010 (file: growth_rates.xlsx)
##
## 		!!!CAUTION!!!: if we use more recent census data (i.e. ADJADM == 1, 
##			see below), it must be the rate from the more recent census data to 
##			2010 (i.e. from 2009 to 2010 for Vietnam) - so calculation, for 
##			VNM, is for one year:
GR_urb_10 = 1.543728
GR_rur_10 = 1.146370

GR_urb_15 = 2.288279
GR_rur_15 = 1.297189


##	Processing flags:

##	Set UNADJUST10 to True if we want to produce a map adjusted to UN totals 
##		for 2010, False otherwise:
UNADJUST10 = True

## If UNADJUST10 == True then we need to provide the UN total population for 
##		2010 - needed if you want to adjust map for U.N. esimates for 
##		2010.  U.N. estimates are from the World Urbanization Prospects 
##		(http://esa.un.org/unpd/wup/index.htm)
UNPOP10 = 158423000

##	Set UNADJUST15 to True if we want to produce a map adjusted to UN totals 
##		for 2015, False otherwise:
UNADJUST15 = True

# UN total population for 2015:
UNPOP15 = 179791000


##	Alternative census folder:
##		This is a census folder other than the one used for the RF model
##		and predict_density.img file.  If not an empty string then the census
##		file in this folder will override the one in the /Census folder!  
##		Usually this folder would start with a "!" as you wouldn't have 
##		wanted it processed.  This is useful if you used a much more detailed, 
##		but older census file to generate RF model predictions and want to 
##		use more concurrent data for population	estimates:
census_folder = ""
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


##	Input projection is configured for each country based on the most
##		appropriate output for distance/area considerations...
##	Final projection for population maps:
final_prj = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"


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

##	Set the output coordinate system and extent to that of our buffer:
arcpy.env.outputCoordinateSystem = final_prj


##	END: Import packages and set GeoProcessing environment
#####



#####
##	BEGIN: Define utility functions

##	Create a directory function to check for and create data directories if they don't exist:
def ensure_dir(d):
	if not os.path.isdir(d):
		os.makedirs(d)
	return d


##	END: Define utility functions
#####



#####
##	BEGIN: Data pre-processing for needed datasets

##	Datasets:
##	Estimated population density from randomForest to be used as density 
##		weightings (gridx) from the old model:
popdensity_weighting = output_path + "predict_density.img"

##	Population data
if census_folder != "":
	census_path = data_path + census_folder + "/"
else:
	census_path = data_path + "Census/"

##	Instead of hardcoding it above we'll just pull in the only shapefile
##		that should be in the directory:
dataset_name = os.path.basename( glob.glob( census_path + "*.shp" )[0] )
adminpop = census_path + dataset_name


##	Get a list of all existing folders in the country's data folder:
dataset_folders = os.walk(data_path).next()[1]


##	NOTE: To properly set the extent we need to query the description of 
##		the shapefile by hand, for some reason setting the extent directly
##		to the name of a dataset doesn't work...
desc = arcpy.Describe(adminpop)
arcpy.env.extent = desc.Extent


##	Set the workspace for data creation and processing:
arcpy.env.workspace = ensure_dir(tmp_path)


##	BEGIN: Data pre-processing for needed datasets
#####



#####
##	BEGIN: Population map creation


print ("PREPROCESS: Finalized density weights...")

##	Project our predicted density weighting file back to GCS1984:
in_path = popdensity_weighting
out_name = "predict_density_prj.tif"

out_path = output_path + "tmp/" + out_name

data_desc = arcpy.Describe(in_path)
input_prj = data_desc.SpatialReference.exportToString()

arcpy.ProjectRaster_management(in_path,out_path,final_prj,"BILINEAR","0.0008333","#","#",input_prj)

##	Set our snapping environment to our new, reprojected and 
##		prediction_density layer:
popdensity_weighting_final = out_path
arcpy.env.snapRaster = out_path



if (arcpy.Exists("admin") == False): arcpy.MakeFeatureLayer_management(adminpop, "admin")
arcpy.FeatureClassToFeatureClass_conversion('admin', arcpy.env.workspace, 'admin_Union')
if (arcpy.Exists("admin_Union") == False): arcpy.MakeFeatureLayer_management("/admin_Union.shp", "admin_Union")

if (len(arcpy.ListFields('admin_Union', "POP")) == 0):
	arcpy.AddField_management('admin_Union', 'POP', 'Double')
arcpy.CalculateField_management('admin_Union', 'POP', ' [ADMINPOP] ', '')



##	We need the landcover file for the country projected into our output projection:
if not ("Landcover" in dataset_folders):
	print("ERROR:  No \"Landcover\" folder found!  This is required and indicates you possibly did not run the \"Data Preparation, R.r\" script or specify configuration options correctly in this Python processing script!")
	exit()
else :
	dataset_folders.remove("Landcover")

##	Clip and project our land cover raster:
dataset_folder = "Landcover"
print("PREPROCESS:  " + dataset_folder)

##	Instead of hardcoding it above we'll just pull in the only TIF or IMG 
##		file that should be in the directory:
dataset_path = (glob.glob( data_path + dataset_folder + "/" + "*.img" ) + glob.glob( data_path + dataset_folder + "/" + "*.tif" ))
if dataset_path:
	dataset_name = os.path.basename( dataset_path[0] )
else:
	print("ERROR:  No land cover file found!")
	exit()

in_path = data_path + dataset_folder + "/" + dataset_name

output_name = "landcover_popmap.tif"

out_path = ensure_dir(data_path + dataset_folder + "/Derived/")
landcover_path = out_path + output_name

data_desc = arcpy.Describe(in_path)
input_prj = data_desc.SpatialReference.exportToString()

if not os.path.isfile(landcover_path) or not skip_existing:
	arcpy.ProjectRaster_management(in_path,landcover_path,final_prj,"NEAREST","0.000833","#","#",input_prj)

landcover = Raster(landcover_path)


##	Population redistribution procedure using population from "admin_Union.shp" 
##		converted to a raster and redistributed according to our weights:
print ("POPMAP: Begin creating population redistribution...")
print ("POPMAP: Using census file - " + adminpop)


gridp = arcpy.FeatureToRaster_conversion('admin_Union', 'POP', 'gridp.img', '0.0008333')
gridy = arcpy.sa.ZonalStatistics('admin_Union', 'ADMINID', popdensity_weighting_final, 'SUM', 'DATA')

popmap = gridp * Raster(popdensity_weighting_final) / gridy
popmap.save(output_path + "/" + country + "_popmap_v" + version + ".tif")


print ("POPMAP: Calculating POPMAP10!")
popmap10 = popmap * (landcover != 190.0) * GR_rur_10 + popmap * (landcover == 190.0) * GR_urb_10
popmap10.save(output_path + country + "_popmap10_v" + version + ".tif")

print ("POPMAP: Calculating POPMAP15!")
popmap15 = popmap * (landcover != 190.0) * GR_rur_15 + popmap * (landcover == 190.0) * GR_urb_15
popmap15.save(output_path + country + "_popmap15_v" + version + ".tif")


if UNADJUST10:
	print ("POPMAP: Calculating POPMAP10Adj!")
	zonsum10 = ZonalStatistics(adminpop,"ISO", popmap10, "SUM","DATA")
	const10 = zonsum10 * 0 + UNPOP10
	popmap10adj = popmap10 * (const10 / zonsum10)
	popmap10adj.save(output_path + country + "_popmap10adj_v" + version + ".tif")

if UNADJUST15:
	print ("POPMAP: Calculating POPMAP15Adj!")
	zonsum15 = ZonalStatistics(adminpop,"ISO", popmap15,"SUM","DATA")
	const15 = zonsum15 * 0 + UNPOP15
	popmap15adj = popmap15 * (const15 / zonsum15)
	popmap15adj.save(output_path + country + "_popmap15adj_v" + version + ".tif")


print ("POPMAP: Completed!")


print ("COMPLETED:  Succesffully created population map outputs!")

##	END: Population map creation
#####
