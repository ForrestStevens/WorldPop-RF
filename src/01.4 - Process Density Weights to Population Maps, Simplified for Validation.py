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
##	Population/census data shapefile should be updated with columns
##		'ADMINID', 'ADMINPOP', 'YEARPOP' and 'ISO'
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

##	Parse main configuration file, which will set the country and root_path
##		variables:
execfile("01.0 - Configuration.py.r")

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

data_path = root_path + "/data/" + country + "/"
output_path = root_path + "/output/" + country + "/"
tmp_path = output_path + "tmp/"


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


##	Get a list of all existing folders in the country's data folder:
dataset_folders = os.walk(data_path).next()[1]


##	Check to see if an alternative census folder exists in the folder list
##		and use it for processing our maps.
##	TODO:  Eventually this needs to be pulled from the Metadata.R file
##		instead of a hard coded directory presence.
if ("! New Census" in dataset_folders):
	census_folder = "! New Census"
else:
	census_folder = ""


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

print("PREPROCESS: Finalized density weights...")

##	Project our predicted density weighting file back to GCS1984:
in_path = popdensity_weighting
out_name = "predict_density_prj.tif"

out_path = output_path + "tmp/" + out_name

data_desc = arcpy.Describe(in_path)
input_prj = data_desc.SpatialReference.exportToString()

arcpy.ProjectRaster_management(in_path,out_path,final_prj,"BILINEAR","0.0008333","#","#",input_prj)

##	Set our snapping environment to our new, reprojected and
##		majority resampled land cover file:
popdensity_weighting_final = out_path
arcpy.env.snapRaster = out_path



##	Create a temporary copy of our census data with an appropriate
##		population field to sum and distribute:
if (arcpy.Exists("admin") == False): arcpy.MakeFeatureLayer_management(adminpop, "admin")
arcpy.FeatureClassToFeatureClass_conversion('admin', arcpy.env.workspace, 'admin_Union')
if (arcpy.Exists("admin_Union") == False): arcpy.MakeFeatureLayer_management("/admin_Union.shp", "admin_Union")

if (len(arcpy.ListFields('admin_Union', "POP")) == 0):
	arcpy.AddField_management('admin_Union', 'POP', 'Double')
arcpy.CalculateField_management('admin_Union', 'POP', ' [ADMINPOP] ', '')



##	Population redistribution procedure using population from "admin_Union.shp"
##		converted to a raster and redistributed according to our weights:
print("POPMAP: Begin creating population redistribution...")
print("PPP: Using census file - " + dataset_name)


##	Read the first row of the attribute table from our census data to
##		determine the year of our census data:
rows = arcpy.SearchCursor(adminpop)
row = rows.next()
census_year = row.YEARPOP


gridp = arcpy.FeatureToRaster_conversion('admin_Union', 'POP', 'gridp.img', '0.0008333')
gridy = arcpy.sa.ZonalStatistics('admin_Union', 'ADMINID', popdensity_weighting_final, 'SUM', 'DATA')

print("PPP: Calculating People Per Pixel " + str(census_year))
popmap = gridp * Raster(popdensity_weighting_final) / gridy

##	NOTE: So there's a bug in the arcpy raster optimization that 
##		often, even though you're specifying the .save() option will not
##		actually write it out to disk until it thinks that you're really
##		removing it from memory.  To force this to happen across versions
##		of ArcGIS and arcpy set the raster object to None like this.
##		Discovered here:
##			http://gis.stackexchange.com/questions/46897/saving-rasters-in-a-python-for-loop-fails-only-on-last-iteration
outPath = output_path + country + "_ppp_v" + rf_version + "_" + str(census_year) + ".tif"
popmap.save(outPath)
popmap = None
popmap = Raster(outPath)


print("POPMAP: Completed!")


print("COMPLETED:  Succesffully created population map outputs!")

##	BEGIN: Population map creation
#####
