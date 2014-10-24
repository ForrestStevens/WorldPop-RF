import arcpy
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")

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



##	Versioning information:
version = "2b"


##	Set the path to the base data folders:
path = 'D:/Documents/Graduate School/Research/Population/Data'

##	Per-country specifics:
##	For final output name
country_name = "Afghanistan"
country_prefix = "AFG"
aggregation_level = "002"


##	Processing flags:
##		NOTE:  All of these were removed for the validation processing...

##	Below here things should be set according to the path and country
##		code above:

##	Define geoprocessing/arcpy environment:
arcpy.env.workspace = path + '/Modeling Work/data'
arcpy.OverwriteOutput = True
arcpy.env.overwriteOutput = True


##	Datasets:
##	Estimated population density from randomForest to be used as density 
##		weightings (gridx) from the old model:
popdensity_weighting = path + "/Modeling Work/output/" + country_prefix + "/predict_density.img"

##	Population data
adminpop = path + "/GIS/Census Data/" + country_name + "/" + "From Andrea/" + country_prefix + "_popdata.shp"	#with columns 'ADMINID', 'ADMINPOP', 'YEARPOP' and 'ISO'

##	Land cover to 100m using a majority resampling, used as a snapping 
##		raster only so replace with a raster common to both modeling 
##		approaches:
refined_landcover_prj = path + "/MDA Land Cover/Refined/" + country_prefix + "/" + country_prefix + "_lc_mosaic_reclass_missing_filled_rurb_8bit_prjGCS84_100m.img"

##	Our final population density weightings (gridx), projected to GCS84:
popdensity_weighting_final = path + "/Modeling Work/output/" + country_prefix + "/predict_density_prjGCS84_final.img"


##	Set the extents for all processing:
arcpy.env.extent = adminpop

##	Set our snapping environment to our new, reprojected and majority resampled
##		land cover file:
arcpy.env.snapRaster = refined_landcover_prj


################## PROCESSING START #############################################

print ("PREPROCESS: Finalized density weights...")


if (arcpy.Exists("admin") == False): arcpy.MakeFeatureLayer_management(adminpop, "admin")
arcpy.FeatureClassToFeatureClass_conversion('admin', arcpy.env.workspace, 'admin_Union')
if (arcpy.Exists("admin_Union") == False): arcpy.MakeFeatureLayer_management("/admin_Union.shp", "admin_Union")

if (len(arcpy.ListFields('admin_Union', "RECNO")) == 0):
	arcpy.AddField_management('admin_Union', 'RECNO', 'Long')
arcpy.CalculateField_management('admin_Union', 'RECNO', ' [FID] + 1 ', '')
if (len(arcpy.ListFields('admin_Union', "POP")) == 0):
	arcpy.AddField_management('admin_Union', 'POP', 'Double')
arcpy.CalculateField_management('admin_Union', 'POP', ' [ADMINPOP] ', '')



##	Population redistribution procedure using population from "admin_Union.shp" 
##		converted to a raster and redistributed according to our weights:
print ("POPMAP: Begin creating population redistribution...")

gridp = arcpy.FeatureToRaster_conversion('admin_Union', 'POP', 'gridp.img', '0.0008333')
gridy = arcpy.sa.ZonalStatistics('admin_Union', 'RECNO', popdensity_weighting_final, 'SUM', 'DATA')

popmap = gridp * Raster(popdensity_weighting_final) / gridy
popmap.save(path + "/Modeling Work/output/" + country_prefix + "/" + country_prefix + "_popmap_admin" + aggregation_level + "_v" + version + ".tif")

print ("POPMAP: Completed!")



print ("COMPLETED:  Succesffully created population map outputs!")
