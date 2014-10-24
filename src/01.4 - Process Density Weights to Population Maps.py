####
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

##	Parse main configuration file, which will set the country and root_path
##		variables:
execfile("01.0 - Configuration.py.r")


##	Round to whole population counts?
round_counts = False



##	TODO: Eventually this can be set in the Metadata.r file and pulled
##		via JSON but for now we will set it here, you must make
##		sure it matches the versioning in the Metadata.r file:


##	Country specific population-specific variables:
##		Growth rates are estimated for urban and rural areas for the years
##		included in GR_years. (file: growth_rates.xlsx)
##
## 		!!!CAUTION!!!: if we use more recent census data (see the
##			census_folder option below), it must be the rate from the more
##			recent census data to	2010 (i.e. from 2009 to 2010 for Vietnam) -
##			so calculation, for	VNM, is for one year:
GR_years = [2010, 2015, 2020]
GR_urb = [0.961943079,1.05992696,1.159165095]
GR_rur = [0.986294783,1.020915731,1.039978459]


##	Processing flags:

##	Set UNADJUST to True if we want to produce a map adjusted to UN totals
##		for 2010, False otherwise:
UNADJUST = [True, True, True]

## If UNADJUST == True then we need to provide the UN total population for
##		that year - needed if you want to adjust map for U.N. esimates.
##		U.N. estimates are from the World Urbanization Prospects
##		(http://esa.un.org/unpd/wup/index.htm)
UNPOP = [5788000, 6213000, 6603000]


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

data_path = root_path + "/data/" + country + "/"
output_path = root_path + "/output/" + country + "/"
tmp_path = output_path + "/tmp/"


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
import time

import arcpy
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *


##	Configure arcpy:

# Set the compression environment to LZ77:
arcpy.env.compression = "LZ77"

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
##		weightings (gridx) from the old model, a file that should be in the
##		directory (we changed formats between version numbers so check for both:
dataset_path = glob.glob( output_path + "predict_density.*" )
if dataset_path:
	popdensity_weighting = dataset_path[0]
else:
	print("ERROR:  No \"predict_density\" TIF or IMG found in the output folder!  You first need to run the 1.3 R script!")
	exit()



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
##		prediction_density layer:
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
print("POPMAP: Begin creating population redistribution...")


print("PPP: Using census file - " + adminpop)

gridp = arcpy.FeatureToRaster_conversion('admin_Union', 'POP', 'gridp.tif', '0.0008333')
gridy = arcpy.sa.ZonalStatistics('admin_Union', 'ADMINID', popdensity_weighting_final, 'SUM', 'DATA')


##	Read the first row of the attribute table from our census data to
##		determine the year of our census data:
rows = arcpy.SearchCursor(adminpop)
row = rows.next()
census_year = row.YEARPOP


##	Calculate population map for census year:
print("PPP: Calculating People Per Pixel " + str(census_year))
if (round_counts):
	##	Int() just truncates to integer values (note any populationcounts >
	##		greater than 2,147,483,647 (maximum size determined by 2^31-1)
	##		will be set to NoData:
	popmap = Int((gridp * Raster(popdensity_weighting_final) / gridy) + 0.5)
else:
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


i = 0
for popyear in GR_years:
	if popyear != census_year:
		print("PPP: Calculating People Per Pixel " + str(popyear))
		if (round_counts):
			##	Int() just truncates to integer values (note any populationcounts >
			##		greater than 2,147,483,647 (maximum size determined by 2^31-1)
			##		will be set to NoData:
			popmap_year = Int((popmap * (landcover != 190.0) * GR_rur[i] + popmap * (landcover == 190.0) * GR_urb[i]) + 0.5)
		else:
			popmap_year = popmap * (landcover != 190.0) * GR_rur[i] + popmap * (landcover == 190.0) * GR_urb[i]

		outPath = output_path + country + "_ppp_v" + rf_version + "_" + str(GR_years[i]) + ".tif"
		popmap_year.save(outPath)
		popmap_year = None
		popmap_year = Raster(outPath)
	else:
		popmap_year = popmap

	if UNADJUST[i]:
		print("PPP: Calculating People Per Pixel, UN Adjusted " + str(popyear))
		zonsum = ZonalStatistics('admin_Union',"ISO", popmap_year, "SUM","DATA")
		const = zonsum * 0 + UNPOP[i]

		if (round_counts):
			##	Int() just truncates to integer values (note any populationcounts >
			##		greater than 2,147,483,647 (maximum size determined by 2^31-1)
			##		will be set to NoData:
			popmap_year_adj = Int( (popmap_year * (const / zonsum)) + 0.5 )
		else:
			popmap_year_adj = popmap_year * (const / zonsum)

		outPath = output_path + country + "_ppp_v" + rf_version + "_" + str(GR_years[i]) + "_UNadj.tif"
		popmap_year_adj.save(outPath)
		popmap_year_adj = None
		popmap_year_adj = Raster(outPath)

	i += 1

print("PPP: Completed!")



##	Now we need to do the same calculations but on the projected data,
##		creating output in people per hectare instead of people per pixel
##		and not converting back to geographic coordinates:

##	Set our snapping environment to the prediction density weighting layer:
popdensity_weighting_final = popdensity_weighting
arcpy.env.snapRaster = popdensity_weighting_final


##	Set the output coordinate system and extent to that of our buffer:
desc = arcpy.Describe(popdensity_weighting_final)
arcpy.env.extent = desc.Extent

final_prj = desc.SpatialReference.exportToString()
arcpy.env.outputCoordinateSystem = final_prj


adminpop_prj = census_path + "Derived/census.shp"
if (arcpy.Exists("admin_prj") == False): arcpy.MakeFeatureLayer_management(adminpop_prj, "admin_prj")
arcpy.FeatureClassToFeatureClass_conversion('admin_prj', arcpy.env.workspace, 'admin_Union_prj')

arcpy.MakeFeatureLayer_management("/admin_Union_prj.shp", "admin_Union_prj")

if (len(arcpy.ListFields('admin_Union_prj', "POP")) == 0):
	arcpy.AddField_management('admin_Union_prj', 'POP', 'Double')
arcpy.CalculateField_management('admin_Union_prj', 'POP', ' [ADMINPOP] ', '')
if (arcpy.Exists("admin_Union_prj") == False): arcpy.MakeFeatureLayer_management("/admin_Union_prj.shp", "admin_Union_prj")


##	Load our already projected landcover file:
out_path = ensure_dir(data_path + dataset_folder + "/Derived/")
landcover_path = out_path + "landcover.tif"


##	Begin processing people per hectare using the same output but assuming
##		our projected data have a linear unit of meters and a square pixel
##		output size of 100m, snapped to the prediction density layer.
print("PPHa: Using census file - " + adminpop)


gridp = arcpy.FeatureToRaster_conversion('admin_Union_prj', 'POP', 'gridp.img', '100')
gridy = arcpy.sa.ZonalStatistics('admin_Union_prj', 'ADMINID', popdensity_weighting_final, 'SUM', 'DATA')


##	Read the first row of the attribute table from our census data to
##		determine the year of our census data:
rows = arcpy.SearchCursor(adminpop)
row = rows.next()
census_year = row.YEARPOP


##	Calculate population map for census year:
print("PPHa: Calculating People Per Hectare " + str(census_year))
if (round_counts):
	##	Int() just truncates to integer values (note any populationcounts >
	##		greater than 2,147,483,647 (maximum size determined by 2^31-1)
	##		will be set to NoData:
	popmap = Int( (gridp * Raster(popdensity_weighting_final) / gridy) + 0.5 )
else:
	popmap = gridp * Raster(popdensity_weighting_final) / gridy

outPath = output_path + country + "_pph_v" + rf_version + "_" + str(census_year) + ".tif"
popmap.save(outPath)
popmap = None
popmap = Raster(outPath)


i = 0
for popyear in GR_years:
	if popyear != census_year:
		print("PPHa: Calculating People Per Hectare " + str(popyear))
		if (round_counts):
			##	Int() just truncates to integer values (note any populationcounts >
			##		greater than 2,147,483,647 (maximum size determined by 2^31-1)
			##		will be set to NoData:
			popmap_year = Int( (popmap * (landcover != 190.0) * GR_rur[i] + popmap * (landcover == 190.0) * GR_urb[i]) + 0.5 )
		else:
			popmap_year = popmap * (landcover != 190.0) * GR_rur[i] + popmap * (landcover == 190.0) * GR_urb[i]

		outPath = output_path + country + "_pph_v" + rf_version + "_" + str(GR_years[i]) + ".tif"
		popmap_year.save(outPath)
		popmap_year = None
		popmap_year = Raster(outPath)
	else:
		popmap_year = popmap

	if UNADJUST[i]:
		print("PPHa: Calculating People Per Hectare, UN Adjusted " + str(popyear))
		zonsum = ZonalStatistics('admin_Union_prj',"ISO", popmap_year, "SUM","DATA")
		const = zonsum * 0 + UNPOP[i]

		if (round_counts):
			##	Int() just truncates to integer values (note any populationcounts >
			##		greater than 2,147,483,647 (maximum size determined by 2^31-1)
			##		will be set to NoData:
			popmap_year_adj = Int( (popmap_year * (const / zonsum)) + 0.5 )
		else:
			popmap_year_adj = popmap_year * (const / zonsum)

		outPath = output_path + country + "_pph_v" + rf_version + "_" + str(GR_years[i]) + "_UNadj.tif"
		popmap_year_adj.save(outPath)
		popmap_year_adj = None
		popmap_year_adj = Raster(outPath)

	i += 1

print("PPHa: Completed!")



print("COMPLETED:  Succesffully created population map outputs!")

##	END: Population map creation
#####