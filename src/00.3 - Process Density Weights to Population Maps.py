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


##	Country specific population-specific variables:
##		GROWTH RATES est for 2010 and 2015 listed here. In this example, listed as 1 bc cntry population census is 2010 and estimated for 2010 (file: growth_rates.xlsx)
## 		!!!CAUTION!!!: if we use more recent census data (i.e. ADJADM == 1, see below), it must be the rate from the more recent census data to 2010
##	(i.e. from 2009 to 2010 for Vietnam) - so calculation, for VNM, is for one year
Expr_urb_10 = 1.038524
Expr_urb_15 = 1.294727

Expr_rur_10 = 1.022551
Expr_rur_15 = 1.172104

##	r = city growth rate: urban growth rate from the year of city pop (from GRUMP shp) to the year of census pop data being used in YEARPOP
r = 1.494363


##	Processing flags:

##	If we want to use city pop estimates from GRUMP: 1, otherwise: 0 - do this based on preprocessing step for cks...if admin units are more spatially refined than GRUMP extents
GRUMPPOP = 0

##	1 if we want to produce a map adjusted to UN totals for 2010, 0 otherwise
UNADJUST10 = 1
# UN total population for 2010 - needed if you want to adjust map for U.N. esimates for 2010 - U.N. estimates are from the World Urbanization Prospects (http://esa.un.org/unpd/wup/index.htm)
UNPOP10 = 31412000

# 1 if we want to produce a map adjusted to UN totals for 2015, 0 otherwise
UNADJUST15 = 1
# UN total population for 2015
UNPOP15 = 36735000


##	For adjustment to coarser resolution but more recent census data (e.g. Vietnam):
##	1 if we want to produce a map adjusted to a coarser admin level, 0 otherwise
ADJADM = 0

if (ADJADM == 1):
	adminpopAdj = "/VNM/VNM-popdata-provinces/VNM-popdata-provinces.shp"	#Additional population data shapefile with columns 'ADMINID', 'ADMINPOP', 'YEARPOP' and 'ISO' that provides the alternative census polygons
	Expr_urb_adj = 1.407057 #GROWTH RATES from old census date to the more recent census (i.e. from 1999 to 2009 for Vietnam) - the first map estimated in process
	Expr_rur_adj = 1.022755



##	Below here things should be set according to the path and country
##		code above:

##	Define geoprocessing/arcpy environment:
arcpy.env.workspace = path + '/Modeling Work/data'
arcpy.OverwriteOutput = True
arcpy.env.overwriteOutput = True


##	Datasets:
##	Land cover data:
##		Note, the first step will be to reproject this data using majority
##		resampling to a cell-size of 0.000833 degrees, GCS_1984 coordinates
##		and then use this as our snapping raster for all further processing...
refined_landcover = path + "/MDA Land Cover/Refined/" + country_prefix + "/" + country_prefix + "_lc_mosaic_reclass_missing_filled_rurb_8bit.img"


##	Estimated population density from randomForest to be used as density 
##		weightings (gridx) from the old model:
popdensity_weighting = path + "/Modeling Work/output/" + country_prefix + "/predict_density.img"

##	Population data
adminpop = path + "/GIS/Census Data/" + country_name + "/" + country_prefix + "_popdata.shp"	#with columns 'ADMINID', 'ADMINPOP', 'YEARPOP' and 'ISO'

####	TODO:  Special case for Afghanistan, we'll just use Andrea's file:
adminpop = path + "/GIS/Census Data/" + country_name + "/From Andrea/" + country_prefix + "_popdata.shp"	#with columns 'ADMINID', 'ADMINPOP', 'YEARPOP' and 'ISO'


##	GRUMP	
GRURB = path + "/GIS/GRUMP/af_as_lac_urban.shp"


##	Set the extents and snapping raster for all processing:
arcpy.env.extent = adminpop



################## PROCESSING START #############################################

##	Project our land cover to 100m using a majority resampling:
refined_landcover_prj = path + "/MDA Land Cover/Refined/" + country_prefix + "/" + country_prefix + "_lc_mosaic_reclass_missing_filled_rurb_8bit_prjGCS84_100m.img"

if (arcpy.Exists(refined_landcover_prj) == False):
	arcpy.ProjectRaster_management(refined_landcover, refined_landcover_prj, "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]","MAJORITY","0.000833","#","#","GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]")

print ("PREPROCESS: Projected land cover...")


##	Set our snapping environment to our new, reprojected and majority resampled
##		land cover file:
arcpy.env.snapRaster = refined_landcover_prj


##	Project our prediction density weightings:
popdensity_weighting_prj = path + "/Modeling Work/output/" + country_prefix + "/predict_density_prjGCS84.img"

if (arcpy.Exists(popdensity_weighting_prj) == False):
	arcpy.ProjectRaster_management(popdensity_weighting, popdensity_weighting_prj,"GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]","NEAREST","0.0008333","#","#","PROJCS['World_Mollweide',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Mollweide'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',0.0],UNIT['Meter',1.0]]")

print ("PREPROCESS: Projected predicted density weights...")

##	Burn in any water from our land cover as 0 in our density weightings:
popdensity_weighting_final = path + "/Modeling Work/output/" + country_prefix + "/predict_density_prjGCS84_final.img"
if (arcpy.Exists(popdensity_weighting_final) == False):
        output = (Raster(refined_landcover_prj) != 210) * Raster(popdensity_weighting_prj)
        output.save(popdensity_weighting_final)

print ("PREPROCESS: Finalized density weights...")


##	TODO: Fix the GRUMP processing to be raster-based...
#if (GRUMPPOP == 1):
#	##	Create shapefile lc_temp (land cover) 
#	arcpy.RasterToPolygon_conversion(refined_landcover_prj, 'lc_temp', 'NO_SIMPLIFY') #CONVERTING RASTER LC TO POLYGON
#	if (arcpy.Exists("lc_temp") == False): arcpy.MakeFeatureLayer_management("/lc_temp.shp", "lc_temp")
#	
#	##	Create file "GRURBsel.shp" (GRUMP urban extents clipped to the country extent) 
#	if (arcpy.Exists("GRURB") == False): arcpy.MakeFeatureLayer_management(GRURB, "GRURB") #USING THE GLOBAL GRUMP FILE AND CUTTING IT DOWN TO COUNTRY SIZE 
#	if (arcpy.Exists("admin") == False): arcpy.MakeFeatureLayer_management(adminpop, "admin")
#	arcpy.SelectLayerByLocation_management('GRURB', 'HAVE_THEIR_CENTER_IN', 'admin', '', 'NEW_SELECTION')
#	arcpy.FeatureClassToFeatureClass_conversion("GRURB", arcpy.env.workspace, 'GRURBsel')
#	##	"SelectLayerByAttribute_management" CREATES, UPDATES, OR REMOVES A SELECTION ON THE LAYER OR TABLE VIEW USING AN ATTRIBUTE	QUERY, IN THIS CASE "CLEARS"
#	arcpy.SelectLayerByAttribute_management('GRURB', 'CLEAR_SELECTION', ' ') 
#
#	#CITYPOP for the right year (adjust from 2000 to the year of census data) # FOR VNM IT IS 2000-1999
#	if (arcpy.Exists("GRURBsel") == False): arcpy.MakeFeatureLayer_management("/GRURBsel.shp", "GRURBsel")
#	arcpy.AddField_management('GRURBsel', 'CITYPOP', 'Double')
#	arcpy.CalculateField_management('GRURBsel', 'CITYPOP', "[POP] * " + str(r), '')
#
#	print ("GRUMP: Successfully created GRUMP polygons that intersect with urban land cover!")
#
#	##	Create file "APURB.shp" (detailed urban extents from the refined land cover map)
#	arcpy.SelectLayerByAttribute_management('lc_temp', 'NEW_SELECTION', ' GRIDCODE = 190') # SELECTING BY GRIDCODE '190' THE URBAN AREA
#
#	arcpy.FeatureClassToFeatureClass_conversion('lc_temp', arcpy.env.workspace, 'APURB')
#	if (arcpy.Exists("APURB") == False): arcpy.MakeFeatureLayer_management("/APURB.shp", "APURB")
#	arcpy.SelectLayerByAttribute_management('lc_temp', 'CLEAR_SELECTION', ' ')
#
#	print ("GRUMP: Successfully intersected and selected GRUMP/urban polgyons!")
#
#	##	Integrate city population numbers from GRUMP urban extents (file "GRURBsel.shp") (C1 to C7 are temporary shapefiles):
#	arcpy.SelectLayerByAttribute_management('GRURBsel', 'NEW_SELECTION', ' CITYPOP <= 1')
#	arcpy.DeleteFeatures_management('GRURBsel')
#	arcpy.SelectLayerByLocation_management('APURB', 'INTERSECT', 'GRURBsel', '', 'NEW_SELECTION')
#
#	arcpy.FeatureClassToFeatureClass_conversion('APURB', arcpy.env.workspace, 'C1')
#	if (arcpy.Exists("C1") == False): arcpy.MakeFeatureLayer_management("/C1.shp", "C1")
#
#	arcpy.SpatialJoin_analysis('C1', 'GRURBsel', 'C2', 'JOIN_ONE_TO_ONE')
#	if (arcpy.Exists("C2") == False): arcpy.MakeFeatureLayer_management("/C2.shp", "C2")
#
#	arcpy.Intersect_analysis('C2;admin', 'C3', 'ALL', '', '')
#	if (arcpy.Exists("C3") == False): arcpy.MakeFeatureLayer_management("/C3.shp", "C3")
#
#	arcpy.Dissolve_management('C3', 'C4', 'ADMINID;ADMINPOP;YEARPOP;ISO;ISOURBID;CITYPOP')
#	if (arcpy.Exists("C4") == False): arcpy.MakeFeatureLayer_management("/C4.shp", "C4")
#
#	arcpy.CalculateAreas_stats('C4', 'C5')
#	if (arcpy.Exists("C5") == False): arcpy.MakeFeatureLayer_management("/C5.shp", "C5")
#
#	arcpy.AddField_management('C5', 'POLYAREA', 'Double')
#	arcpy.CalculateField_management('C5', 'POLYAREA', ' [F_AREA] * 1000000 ', '')
#	arcpy.Statistics_analysis('C5', 'st1.dbf', 'POLYAREA sum', 'ISOURBID')
#	arcpy.AddJoin_management('C5', 'ISOURBID', 'st1.dbf', 'ISOURBID', 'KEEP_ALL')
#
#	arcpy.FeatureClassToFeatureClass_conversion('C5', arcpy.env.workspace, 'C6')
#	if (arcpy.Exists("C6") == False): arcpy.MakeFeatureLayer_management("/C6.shp", "C6")
#
#	arcpy.AddField_management('C6', 'POLYPOP', 'Double')
#	arcpy.CalculateField_management('C6', 'POLYPOP', ' ( [POLYAREA] / [st1_sum_PO] ) * [CITYPOP] ', '')
#	arcpy.DeleteField_management('C6', 'st1_OID;st1_ISOURB;st1_FREQUE')
#	arcpy.Statistics_analysis('C6', 'st2.dbf', 'POLYPOP sum', 'ADMINID')
#	arcpy.AddJoin_management('C6', 'ADMINID', 'st2.dbf', 'ADMINID', 'KEEP_ALL')
#
#	arcpy.FeatureClassToFeatureClass_conversion('C6', arcpy.env.workspace, 'C7')
#	if (arcpy.Exists("C7") == False): arcpy.MakeFeatureLayer_management("/C7.shp", "C7")
#
#	arcpy.AddField_management('C7', 'URBPOP', 'Double')
#	arcpy.CalculateField_management('C7', 'URBPOP', ' [st2_sum_PO] ', '')
#	arcpy.DeleteField_management('C7', 'st2_OID;st2_ADMINI;st2_FREQUE;st2_sum_PO')
#	arcpy.AddJoin_management('admin', 'ADMINID', 'C7', 'ADMINID', 'KEEP_ALL')
#
#	arcpy.FeatureClassToFeatureClass_conversion('admin', arcpy.env.workspace, 'admin2')
#	if (arcpy.Exists("admin2") == False): arcpy.MakeFeatureLayer_management("/admin2.shp", "admin2")
#
#	arcpy.AddField_management('admin2', 'RURPOP', 'Double')
#	arcpy.CalculateField_management('admin2', 'RURPOP', ' [ADMINPOP] - [C7_URBPOP] ', '')
#
#	print ("GRUMP: Successfully aggregated urban and rural population numbers into totals!")
#
#	##	Create "admin_Union.shp" with column "POP"
#	arcpy.Union_analysis('admin2;C7', 'admin_Union')
#	if (arcpy.Exists("admin_Union") == False): arcpy.MakeFeatureLayer_management("/admin_Union.shp", "admin_Union")
#
#	if (len(arcpy.ListFields('admin_Union', "RECNO")) == 0):
#		arcpy.AddField_management('admin_Union', 'RECNO', 'Long')
#	arcpy.CalculateField_management('admin_Union', 'RECNO', ' [FID] + 1 ', '')
#
#	if (len(arcpy.ListFields('admin_Union', "POP")) == 0):
#		arcpy.AddField_management('admin_Union', 'POP', 'Double')
#	
#	arcpy.SelectLayerByAttribute_management('admin_Union', 'NEW_SELECTION', ' CITYPOP = 0')
#	arcpy.CalculateField_management('admin_Union', 'POP', ' [RURPOP] ', '')
#
#	arcpy.SelectLayerByAttribute_management('admin_Union', 'NEW_SELECTION', ' CITYPOP > 0')
#	arcpy.CalculateField_management('admin_Union', 'POP', ' [POLYPOP] ', '')
#
#	arcpy.SelectLayerByAttribute_management('admin_Union', 'NEW_SELECTION', ' POP < 0')
#	arcpy.CalculateField_management('admin_Union', 'POP', ' 0 ', '')
#
#	arcpy.SelectLayerByAttribute_management('admin_Union', 'CLEAR_SELECTION', ' ')
#
#	print ("GRUMP: Completed!")


if (GRUMPPOP == 0):
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

if (ADJADM == 0):
	popmap = gridp * Raster(popdensity_weighting_final) / gridy
	popmap.save(path + "/Modeling Work/output/" + country_prefix + "/" + country_prefix + "_popmap_v" + version + ".tif")

if (ADJADM == 1):
	popmaptemp = gridp * Raster(popdensity_weighting_final) / gridy
	popmaptemp.save("popmaptemp.img")

print ("POPMAP: Completed!")



##	If ADJADM = 1 this part of the code is run to adjust population values and re-estimate for final product
if (ADJADM == 1):
	print ("ADJUSTED: Begin adjusting population map for coarser/recent census data...")
       
	##	Project forward to the year of the more recent census data (e.g. 1999 to 2009 for Vietnam):
	popmaptempAdj = popmaptemp * (Raster(refined_landcover_prj) != 190) * Expr_rur_adj + popmaptemp * (Raster(refined_landcover_prj) == 190) * Expr_urb_adj
	popmaptempAdj.save("popmaptempAdj.img")


	##	Adjustment to the population values:
	if (arcpy.Exists("adminpopAdj") == False): arcpy.MakeFeatureLayer_management(adminpopAdj, "adminpopAdj") 
	outZStat = arcpy.sa.ZonalStatistics(adminpopAdj, "ADMINID", popmaptempAdj, "SUM", "DATA")
	
	##	Create a raster from our new/coarse census data:
	gridpAdj = arcpy.FeatureToRaster_conversion(adminpopAdj, 'ADMINPOP', 'gridpAdj.img', '0.0008333')

	##	Now scale our redistributed population by the ratio of coarse/new census data
	##		to the sum across the block of our redistributed data:
	popmap = (gridpAdj / outZStat) * popmaptemp
	popmap.save(path + "/Modeling Work/output/" + country_prefix + "/" + country_prefix + "_popmap_v" + version + ".tif")

	print ("ADJUSTED: Population correctly adjusted for more recent census data:")

##	Project forward to 2010:
popmap10 = popmap * (Raster(refined_landcover_prj) != 190) * Expr_rur_10 + popmap * (Raster(refined_landcover_prj) == 190) * Expr_urb_10
popmap10.save(path + "/Modeling Work/output/" + country_prefix + "/" + country_prefix + "_popmap10_v" + version + ".tif")

print ("PROJECTED: Population projected forward to 2010...")


##	Project forward to 2015:
popmap15 = popmap * (Raster(refined_landcover_prj) != 190) * Expr_rur_15 + popmap * (Raster(refined_landcover_prj) == 190) * Expr_urb_15
popmap15.save(path + "/Modeling Work/output/" + country_prefix + "/" + country_prefix + "_popmap15_v" + version + ".tif")

print ("PROJECTED: Population projected forward to 2015...")


if (UNADJUST10 == 1):
	##	Run zonal stats to create a raster with total population in each pixel
	##		based on population map we estimmated (ISO used as field because 
	##		standardized across all admin units):
	outZStat = arcpy.sa.ZonalStatistics(adminpop, "ISO", popmap10, "SUM", "DATA")
	
	##	Create raster with UN POP10 with total pop/pixel from the U.N.:
	outConstRast = arcpy.sa.CreateConstantRaster(UNPOP10, "INTEGER", '0.0008333', popmap10.extent)

	##	Finish adjustment to U.N. population totals for 2010:
	popmap10adj = popmap10 * (outConstRast / outZStat)
	popmap10adj.save(path + "/Modeling Work/output/" + country_prefix + "/" + country_prefix + "_popmap10_UN_v" + version + ".tif")

	print ("UN ADJUST: Population estimates adjusted to match 2010 UN...")

if (UNADJUST15 == 1): # same process but for 2015
	outZStat = arcpy.sa.ZonalStatistics(adminpop, "ISO", popmap15, "SUM", "DATA")
	
	##	Create raster with UN POP15 with total pop/pixel from the U.N.:
	outConstRast = arcpy.sa.CreateConstantRaster(UNPOP15, "INTEGER", '0.0008333', popmap15.extent)

	##	Finish adjustment to U.N. population totals for 2015:
	popmap15adj = popmap15 * (outConstRast / outZStat)
	popmap15adj.save(path + "/Modeling Work/output/" + country_prefix + "/" + country_prefix + "_popmap15_UN_v" + version + ".tif")

	print ("UN ADJUST: Population estimates adjusted to match 2015 UN...")



print ("COMPLETED:  Succesffully created population map outputs!")
