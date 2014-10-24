require(rgdal)
require(raster)
require(RPyGeo)
require(randomForest)
require(foreign)
##	NOTE: The tcltk library is used to create a progress bar, but the
##		raster package has a progress bar mechanism built in in recent versions.
##		Look up the pbCreate(), pbStep(), and pbClose() functions in the raster
##		package...
require(tcltk)


##	Parameters and defaults:

country_name <- "Afghanistan"
country_prefix <- "AFG"

##	The name of the field containing the census counts to estimate and 
##		dasymetrically redistribute:
popfield <- "EST0910"

#path = "D:/Documents/Graduate School/Research/Population/Data"
path = "C:/Users/forrest/Research/Population/Data"

proj4str_gcs_wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj4str_moll_wgs84 <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"


##	Functions:
my_zonalstats <- function(raster_file=paste(path, "/GIS/DIVA-GIS/", country_name, "/Neighbor Merge/", country_prefix, "_gaz_merge_pop_points_count.img", sep=""), census_file = paste(path, "/GIS/Census Data/", country_name, "/", country_prefix,"_popdata_prj_areas.shp", sep=""), FUN="SUM") {
	##	We could do this all in R, using a zonal raster generated from our census_data polygons, however, the raster package's processing speed is just too slow...  :(  Therefore, we're going to have to rely on ArcGIS to do the processing for us.

	###	First we create a zonal raster based on OBJECTID field:
	#setwd(paste(path, "GIS/DIVA-GIS", country_name, "Neighbor Merge", sep="/"))
	#
	#if (file.exists(paste(country_prefix, "_zones.img", sep=""))) {
	#	zonal_raster <- brick(paste(country_prefix, "_zones.img", sep=""))
	#} else {
	#	zonal_raster <- raster(pop_points_count)
	#	zonal_raster <- rasterize(census_data, zonal_raster, field="OBJECTID", datatype="INT2U", format="HFA", filename=paste(country_prefix, "_zones.img", sep=""), progress="text")
	#	#zonal_raster <- rasterize(census_data, zonal_raster, field="OBJECTID", datatype="INT2U", progress="text")
	#}
	#
	#
	###	Next we'll calculate sum of population points:
	#output <- zonal(x=pop_points_count, zone=zonal_raster, stat="sum", na.rm=TRUE, progress="text")

	##	Instead of the R zonal statistical technique, I'm going to use RPyGeos to
	##		access the Python ArcGIS geoprocessing environments:
	workspace_dir <- "C:/tmp"

	rpygeo_env <- rpygeo.build.env(extensions="Spatial", python.path="C:/Python26/ArcGIS10.0", python.command="python.exe", workspace=workspace_dir, overwriteoutput=1)

	###	Could use arcpy...  Note the added four spaces so that the second and fourth lines of the Python code line up with the indentation inside the try/catch block in the generated Python file.  Check the Python file for syntax errors if you're having problems...
	#rpygeo.geoprocessor(
	#	paste(
	#"import arcpy
	#    arcpy.CheckOutExtension(\"Spatial\")
	#
	#    arcpy.sa.ZonalStatisticsAsTable(\"", census_file, "\",\"OBJECTID\",\"", ,"\",\"", path, "/Modeling Work/output/census_data_pop_points_count_sum.dbf\",\"DATA\",\"SUM\")", 
	#	sep=""),
	#env=rpygeo_env, add.gp=FALSE, clean.up=FALSE, working.directory="C:/tmp")

	##	The overwriteoutput=1 flag doesn't seem to work (the geoprocessor says that the workspace is read-only so we'll just check to see if our output DBF file is there and remove it if so:
	#unlink( paste(workspace_dir, "/my_zonalstats_output.dbf", sep="") )

	##	Or just use the gp envronment...
	rpygeo.geoprocessor(
		paste(
	"gp.ZonalStatisticsAsTable_sa(\"", census_file, "\",\"OBJECTID\",\"", raster_file,"\",\"", workspace_dir, "/my_zonalstats_output.dbf\",\"DATA\",\"", FUN, "\")", 
		sep=""),
		env=rpygeo_env, 
		add.gp=FALSE, 
		clean.up=FALSE, 
		working.directory="C:/tmp"
	)

	##	Check to see if there's an error and bail:
	if (file.exists(paste(workspace_dir, "/rpygeo.msg", sep=""))) {
		print("ERROR:  Error in geoprocessing the zonal statistics... Check the .msg file in the geoprocessing working directory!")
		flush.console()
		return(NULL)
	}

	##	Now load up the DBF and join it to our census data and add the sum to
	##		out census_data with covariates data.frame:
	output_dbf <- foreign::read.dbf("C:/tmp/my_zonalstats_output.dbf")

	return(output_dbf)
}


##	Read in our raster-based covariates:
setwd(paste(path, "GIS/DIVA-GIS", country_name, "Neighbor Merge", sep="/"))

pop_points_count <- brick(paste(country_prefix, "_gaz_merge_pop_points_count.img", sep=""))
pop_points_dist <- brick(paste(country_prefix, "_gaz_merge_pop_points_dist.img", sep=""))
roads_dist <- brick(paste(country_prefix, "_roads_merge_dist.img", sep=""))
roads_prim_dist <- brick(paste(country_prefix, "_roads_merge_prim_dist.img", sep=""))
water_areas_dist <- brick(paste(country_prefix, "_water_areas_dcw_merge_dist.img", sep=""))
water_areas_perm_dist <- brick(paste(country_prefix, "_water_areas_dcw_merge_perm_dist.img", sep=""))
water_lines_dist <- brick(paste(country_prefix, "_water_lines_dcw_merge_dist.img", sep=""))
water_lines_perm_dist <- brick(paste(country_prefix, "_water_lines_dcw_merge_perm_dist.img", sep=""))


##	Load the updated and recoded climate zone data, projected into World Mollweide:
setwd(paste(path, "GIS/Peel Climate Classification/Recoded", sep="/"))
climate_zones <- brick("AFG_koppen_recode_prj.img")
#climate_zones_prj <- projectRaster(from=climate_zones, crs=proj4str_moll_wgs84, method="ngb")


##	Load the HydroSHEDS DEM and derived slope rasters:
setwd(paste(path, "HydroSHEDS/Void-Filled DEM", sep="/"))
dem_elev <- brick(paste(country_prefix, "_dem_prj.img", sep=""))
dem_slope <- brick(paste(country_prefix, "_dem_prj_slope.img", sep=""))


##	Load the individual raster layers corresponding to land cover class masks, aggregated via SUM to 100m and then SUM from a circular moving window of radius 9 pixels:
##	NOTE: We are using the _agg.img not the _focal.img to look at 100m pixel scale.  Recall that the 100m pixel scale _agg.img were aggregated using SUM function from the 25m data, so each 100m pixel consists of 16 individual 25m pixels, and proportion makes up the proportion within that 100m cell. 
setwd(paste(path, "MDA Land Cover/Refined/", country_prefix, sep="/"))
#land_cover_v011_f <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v011_agg_focal.img", sep=""))
#land_cover_v040_f <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v040_agg_focal.img", sep=""))
#land_cover_v130_f <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v130_agg_focal.img", sep=""))
#land_cover_v140_f <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v140_agg_focal.img", sep=""))
#land_cover_v160_f <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v160_agg_focal.img", sep=""))
land_cover_v190_f <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v190_agg_focal.img", sep=""))
#land_cover_v200_f <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v200_agg_focal.img", sep=""))
#land_cover_v210_f <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v210_agg_focal.img", sep=""))
land_cover_v240_f<- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v240_agg_focal.img", sep=""))
land_cover_vSum <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_vNNN_agg_focal_sum.img", sep=""))

land_cover_v011 <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v011_agg.img", sep=""))
land_cover_v040 <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v040_agg.img", sep=""))
land_cover_v130 <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v130_agg.img", sep=""))
land_cover_v140 <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v140_agg.img", sep=""))
land_cover_v160 <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v160_agg.img", sep=""))
land_cover_v190 <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v190_agg.img", sep=""))
land_cover_v200 <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v200_agg.img", sep=""))
land_cover_v210 <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v210_agg.img", sep=""))
land_cover_v240 <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v240_agg.img", sep=""))
land_cover_v190v240 <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v190_v240_agg.img", sep=""))

##	Load the distance-to land cover mask data:
land_cover_v011_dist <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v011_agg_mask_dist.img", sep=""))
land_cover_v190_dist <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v190_agg_mask_dist.img", sep=""))
land_cover_v240_dist <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v240_agg_mask_dist.img", sep=""))
land_cover_v190v240_dist <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v190_v240_agg_mask_dist.img", sep=""))



##	Load the land cover classification, projected into World Mollweide:
setwd(paste(path, "MDA Land Cover/Refined/", country_prefix, sep="/"))
##	readGDAL() isn't an option because the data file is too large... Therefore
##		we need to use the raster package functionality...
#land_cover <- readGDAL(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj.img", sep=""))
land_cover <- brick(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj.img", sep=""))
#land_cover <- raster(paste(country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj.img", sep=""))
#land_cover_prj <- projectRaster(from=land_cover, crs=proj4str_moll_wgs84, method="ngb")


#####
##	If our covariate data already exist load it, if not create it:
setwd(paste(path, "GIS/Census Data", country_name, sep="/"))

if (file.exists(paste(country_prefix, "_popdata_prj_areas_covariates.shp", sep=""))) {
	census_data <- readOGR(".", paste(country_prefix, "_popdata_prj_areas_covariates", sep=""))
} else {

	##	Define our census file for the zonal statistics calculation:
	census_file <- paste(path, "/GIS/Census Data/", country_name, "/", country_prefix,"_popdata_prj_areas.shp", sep="")

	##	Read in our census data with final population counts and area field, F_AREA (in m^2):
	census_data <- readOGR(".", paste(country_prefix, "_popdata_prj_areas", sep=""))

	##	Check to see if we are in Mollweide projection and projects if not:
	#if (projection(census_data) != proj4str_moll_wgs84) {
	#	census_data <- spTransform(census_data, CRS(proj4str_moll_wgs84))
	#}

	##	If you haven't used the "Calculate Areas" tool in ArcGIS to add the F_AREA field
	##		to the shapefile, you theoretically could get them from the readOGR() input
	##		sp data object:
	#areas <- unlist( sapply(slot(census_data, "polygons"), function(x) sapply(slot(x,"Polygons"), slot, "area")) )

	##	Recalculate our F_AREA into hectares:
	census_data$F_AREA_HA <- census_data$F_AREA / 10000

	##	Finally calculate our population density in people per hectare for use as a
	##		dependent variable:
	census_data$POPD_PPHA <- census_data[[popfield]] / census_data$F_AREA_HA

##	NOTE: I could do it this way, using the raster package's functionality, but instead I ended up exporting rasters that are binary masks for each land cover type that are much more efficient to process using the zonal stats functions from ArcGIS' geoprocessor.  See the updated version of this below:

#	##	For each polygon in our dataset we are going to calculate the proportion of each
#	##		covered in our land cover classes to use as predictors:
#
#	##	Create time stamp:
#	start_time = proc.time()[3]
#
#	##	For the progress bar:
#	min_index <- 1
#	max_index <- length(census_data)
#	pb <- tkProgressBar(title = "Generating Land Cover Proportions:", min = min_index, max = max_index, width = 300)
#
#	for (i in 1:length(census_data)) {
#		##	Extract the data for our current polygon:
#		##		NOTE: This is the fastest way I've found, get the extent of the current
#		##			census polygon and create an in-memory RasterLayer object from our 
#		##			land_cover data by using the crop() function.  Then we extract the values
#		##			of that in-memory layer:
#		land_cover_subset <- crop(land_cover, extent(census_data[i,]))
#		land_cover_subset_data <- unlist(extract(land_cover_subset, census_data[i,]))
#
#		##	In our case, 0 values are NAs, so we convert them:
#		land_cover_subset_data[land_cover_subset_data == 0] <- NA
#
#		##	If our na.rm flag is TRUE remove any NAs from the value list:
#		land_cover_subset_data <- land_cover_subset_data[!is.na(land_cover_subset_data)]
#
#		##	Generate counts of each land cover type:
#		counts <- table(land_cover_subset_data)
#		
#		##	Generate a total sum of cells:
#		total <- sum(counts)
#
#		##	Generate the proportion of each land cover within the census block
#		##		and store it in the data.frame:
#		for (j in c("11", "40", "130", "140", "160", "190", "200", "210", "240")) {
#			if(!is.na(counts[j])) {
#				census_data@data[i, paste("PROP_LC_", j, sep="")] <- counts[j] / total
#			} else {
#				census_data@data[i, paste("PROP_LC_", j, sep="")] <- 0
#			}
#		}
#		
#		setTkProgressBar(pb, i, label=paste( round((i - min_index)/(max_index - min_index)*100, 0),"% done"))
#	}
#	close(pb)
#
#	print(paste("Elapsed Time:", proc.time()[3] - start_time, "seconds"))
#
#	##	Confirm that our added fields sum to one:
#	#sum(census_data@data[399,22:30])

	##	For each land cover mask we need to calculate zonal statistics and then calculate the proportion of the cells within the polygon made up by that land cover and add it to our dataset:
	
	for (j in c("011", "040", "130", "140", "160", "190", "200", "210", "240")) {
		raster_file <- paste(path, "/MDA Land Cover/Refined/", country_prefix, "/", country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v", j, "_agg.img", sep="")
		output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="SUM")

		##	NOTE: We need to multiply the count by 16 because these rasters are aggregated to 100m from 25m land cover by a SUM function:
		census_data@data[, paste("PLC_", j, sep="")] <- output_dbf$SUM / (16 * output_dbf$COUNT)
	}


	##	Now we need to aggregate our other covariate data by census block:

	##	Population points count, summed:
	raster_file <- paste(path, "/GIS/DIVA-GIS/", country_name, "/Neighbor Merge/", country_prefix, "_gaz_merge_pop_points_count.img", sep="")
	output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="SUM")
	census_data@data$pp_sum <- output_dbf[,4]

	##	Population points distance, mean:
	raster_file <- paste(path, "/GIS/DIVA-GIS/", country_name, "/Neighbor Merge/", country_prefix, "_gaz_merge_pop_points_dist.img", sep="")
	output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="MEAN")
	census_data@data$pp_dist <- output_dbf[,4]

	##	Roads distance, mean:
	raster_file <- paste(path, "/GIS/DIVA-GIS/", country_name, "/Neighbor Merge/", country_prefix, "_roads_merge_dist.img", sep="")
	output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="MEAN")
	census_data@data$rd_dist <- output_dbf[,4]

	##	Primary roads distance, mean:
	raster_file <- paste(path, "/GIS/DIVA-GIS/", country_name, "/Neighbor Merge/", country_prefix, "_roads_merge_prim_dist.img", sep="")
	output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="MEAN")
	census_data@data$rdp_dist <- output_dbf[,4]

	##	Water areas distance, mean:
	raster_file <- paste(path, "/GIS/DIVA-GIS/", country_name, "/Neighbor Merge/", country_prefix, "_water_areas_dcw_merge_dist.img", sep="")
	output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="MEAN")
	census_data@data$wta_dist <- output_dbf[,4]

	##	Permanent water areas distance, mean:
	raster_file <- paste(path, "/GIS/DIVA-GIS/", country_name, "/Neighbor Merge/", country_prefix, "_water_areas_dcw_merge_perm_dist.img", sep="")
	output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="MEAN")
	census_data@data$wtap_dist <- output_dbf[,4]

	##	Streams/rivers distance, mean:
	raster_file <- paste(path, "/GIS/DIVA-GIS/", country_name, "/Neighbor Merge/", country_prefix, "_water_lines_dcw_merge_dist.img", sep="")
	output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="MEAN")
	census_data@data$wtl_dist <- output_dbf[,4]

	##	Permanent streams/rivers distance, mean:
	raster_file <- paste(path, "/GIS/DIVA-GIS/", country_name, "/Neighbor Merge/", country_prefix, "_water_lines_dcw_merge_perm_dist.img", sep="")
	output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="MEAN")
	census_data@data$wtlp_dist <- output_dbf[,4]

	##	DEM, elevation:
	raster_file <- paste(path, "/HydroSHEDS/Void-Filled DEM/", country_prefix, "_dem_prj.img", sep="")
	output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="MEAN")
	census_data@data$dem_elev <- output_dbf[,4]
	
	##	DEM, slope:
	raster_file <- paste(path, "/HydroSHEDS/Void-Filled DEM/", country_prefix, "_dem_prj_slope.img", sep="")
	output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="MEAN")
	census_data@data$dem_slope <- output_dbf[,4]
	
	##	Peel, Koeppen-Geiger Climate CLassification, majority:
	raster_file <- paste(path, "/GIS/Peel Climate Classification/Recoded/", country_prefix, "_koppen_recode_prj.img", sep="")
	output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="MAJORITY")
	census_data@data$clim_zone <- output_dbf[,4]
	
	
	##	For selected land covers (agriculture, urban, and rural) we need to average the distance-to a 1/NULL mask we generated in ArcGIS (see the "MDA Land Cover/NOTES.txt" file).  We need to calculate zonal statistics:
	for (j in c("011", "190", "240")) {
		raster_file <- paste(path, "/MDA Land Cover/Refined/", country_prefix, "/", country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v", j, "_agg_mask_dist.img", sep="")
		output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="MEAN")
	census_data@data[, paste("DLC_", j, sep="")] <- output_dbf[,4]
	}

	##	Generate proportion and distance to "built" land cover from the original GeoCover data:
	j <- "190_v240"
	raster_file <- paste(path, "/MDA Land Cover/Refined/", country_prefix, "/", country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v", j, "_agg.img", sep="")
	output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="SUM")
	##	NOTE: We need to multiply the count by 16 because these rasters are aggregated to 100m from 25m land cover by a SUM function:
	census_data@data[, paste("PLC_Built", sep="")] <- output_dbf$SUM / (16 * output_dbf$COUNT)
	
	raster_file <- paste(path, "/MDA Land Cover/Refined/", country_prefix, "/", country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v", j, "_agg_mask_dist.img", sep="")
	output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="MEAN")
	census_data@data[, paste("DLC_Built", sep="")] <- output_dbf[,4]
	
	
	##	For selected land covers, we are going to use a focal analysis to 
	##		come up with the proportion of that land cover within ~1km (9
	##		pixel radius) of that particular point:
	
	for (j in c("190", "240")) {
		raster_file <- paste(path, "/MDA Land Cover/Refined/", country_prefix, "/", country_prefix, "_lc_mosaic_reclass_missing_filled_rurb_8bit_prj_v", j, "_agg_focal.img", sep="")
		output_dbf <- my_zonalstats(census_file=census_file, raster_file=raster_file, FUN="SUM")

		##	NOTE: We need to multiply the count by 4048 because these rasters are aggregated to a 9 pixel radius from 100m land cover, previously aggregated from 25m by SUM functions:
		census_data@data[, paste("PLC1k_", j, sep="")] <- output_dbf$SUM / (4048 * output_dbf$COUNT)
	}




	##	Write out our shapefile with covariates:
	writeOGR(census_data, ".", paste(country_prefix, "_popdata_prj_areas_covariates", sep=""), driver="ESRI Shapefile", check_exists=TRUE, overwrite_layer=TRUE)
}

##	Load in our projected census data zonal mask, used as a template for predictions:
census_data_zones <- brick(paste(country_prefix, "_popdata_prj_zones.img", sep=""))
#####

#####
##	Set up response and covariate dataframes for the random forest modeling:
y_data <- census_data@data$POPD_PPHA

##	Full dataset:
#x_data <- census_data@data[,22:48]

### Subset dataset with only primary water and road areas:
#x_data <- census_data@data[,c(22:32, 34, 36, 38:44)]

### Subset dataset without populated places, and just important
###	drivers:
#x_data <- census_data@data[,c(22:30, 34, 36, 38:44)]

### Subset dataset without populated places, and just primary
###	drivers excluding DEM information:
#x_data <- census_data@data[,c(22:30, 34, 36, 38, 41:44)]

### Subset dataset with just populated places sum, and just primary
###	drivers excluding DEM information:
#x_data <- census_data@data[,c(22:31, 34, 36, 38, 41:44)]

### Subset dataset without populated places, and just primary
###	drivers excluding DEM information, and using the "Built" class
###	instead of the v190 and v240 classes separately:
#x_data <- census_data@data[,c(22:26, 28:29, 34, 36, 38, 41:42, 45:46)]

### Subset dataset without populated places, and just primary
###	drivers excluding DEM information, and using the "Built" class
###	instead of the v190 and v240 classes separately, but with a
###	proportion of "urban" and "rural" within 1k of pixel also
###	included:
#x_data <- census_data@data[,c(22:26, 28:29, 34, 36, 38, 41:42, 45:47)]

## Subset dataset without populated places, and just primary
##	drivers excluding DEM information, and using the "Built" class
##	instead of the v190 and v240 classes separately, but with a
##	proportion of "urban" and "rural" within 1k of pixel also
##	included, with PLC_140 and wtap_dist removed (had a negative Inc. %MSE):
x_data <- census_data@data[,c(22:24, 26, 28:29, 34, 38, 41:42, 45:47)]

##	Convert to factors where appropriate:
x_data$clim_zone <- factor(x_data$clim_zone, levels=1:7, labels=as.character(1:7))

###	Subset x_data to land cover alone:
#x_data <- census_data@data[,c(22:30)]

#index <- y_data < 10 	##	Subset data to remove outliers...
index <- y_data > -1 		##	Don't subset data...
y_data <- y_data[index]
x_data <- x_data[index,]



##########
##	Random Forest regression of population density:
set.seed(2002)

##	Check our memory footprint just to make sure:
memory.size()


##	Now we are going tune our randomForest population density regression: 
tuneRF(x=x_data, y=y_data, plot=TRUE, mtryStart=9, ntreeTry=5000, improve=0.0001, stepFactor=1.25, trace=TRUE) 

##mtry = 9  OOB error = 13.94321 
##Searching left ...
##mtry = 8        OOB error = 14.36247 
##-0.03006929 1e-04 
##Searching right ...
##mtry = 11       OOB error = 14.2606 
##-0.02276278 1e-04 
##   mtry OOBError
##8     8 14.36247
##9     9 13.94321
##11   11 14.26060

#	Based on the above I'll use an mtry value of 9:
start_time = proc.time()[3]
popfit = randomForest(x=x_data, y=y_data, mtry=9, nodesize=5, ntree = 5000, na.action = na.omit, importance = TRUE, keep.forest=TRUE)
print(paste("Elapsed Fitting Time:", proc.time()[3] - start_time, "seconds"))
popfit

##Call:
## randomForest(x = x_data, y = y_data, ntree = 5000, mtry = 9,      nodesize = 5, importance = TRUE, keep.forest = TRUE, na.action = na.omit) 
##               Type of random forest: regression
##                     Number of trees: 5000
##No. of variables tried at each split: 9
##
##          Mean of squared residuals: 14.21077
##                    % Var explained: 77.93
					
#summary(popfit)
#plot(popfit)
#
#
###	Check our diagnostics:
###	Recall that running predict on the randomForest object will return OOB predictions:
#predict(popfit)
#importance(popfit)
#
###	Get the variable names of the top 20 predictors:
#names(importance(popfit)[,"%IncMSE"][order(importance(popfit)[,"%IncMSE"], decreasing=TRUE)])[1:20]
#
###	Sort all covariates by they $IncMSE:
#importance(popfit)[order(importance(popfit)[,1], decreasing=TRUE),]
#
#varImpPlot(popfit)
#varUsed(popfit)
#summary(treesize(popfit))
#partialPlot(x=popfit, pred.data=x_data, x.var="PLC_011")
#partialPlot(x=popfit, pred.data=x_data, x.var="DLC_011")
##partialPlot(x=popfit, pred.data=x_data, x.var="PLC_190")
##partialPlot(x=popfit, pred.data=x_data, x.var="DLC_190")
##partialPlot(x=popfit, pred.data=x_data, x.var="PLC_240")
##partialPlot(x=popfit, pred.data=x_data, x.var="DLC_240")
#partialPlot(x=popfit, pred.data=x_data, x.var="PLC_Built")
#partialPlot(x=popfit, pred.data=x_data, x.var="DLC_Built")
#partialPlot(x=popfit, pred.data=x_data, x.var="PLC1k_190")
#partialPlot(x=popfit, pred.data=x_data, x.var="rd_dist")
#partialPlot(x=popfit, pred.data=x_data, x.var="rdp_dist")
#partialPlot(x=popfit, pred.data=x_data, x.var="wtlp_dist")
#partialPlot(x=popfit, pred.data=x_data, x.var="wtap_dist")
#partialPlot(x=popfit, pred.data=x_data, x.var="dem_elev")
#partialPlot(x=popfit, pred.data=x_data, x.var="dem_slope")
#partialPlot(x=popfit, pred.data=x_data, x.var="pp_sum")
#partialPlot(x=popfit, pred.data=x_data, x.var="pp_dist")
#partialPlot(x=popfit, pred.data=x_data, x.var="clim_zone")
#plot(popfit)
#
#
###	For continuous regression, plot observed vs. predicted:
#plot(x=y_data, y=predict(popfit), ylim=c(0,5), xlim=c(0,5))
#plot(x=y_data, y=predict(popfit), ylim=c(0,120), xlim=c(0,120))
#abline(a=0, b=1, lty=2)
#
###	For continuous regression, plot residuals vs. observed:
#plot(x=y_data, y=(y_data - predict(popfit)), xlim=c(0,20))
#abline(a=0, b=0, lty=2)
#
###	For continuous regression, plot residuals vs. fitted:
#plot(x=predict(popfit), y=(y_data - predict(popfit)), xlim=c(0,20))
#abline(a=0, b=0, lty=2)
##########



##########
##	yaImpute imputation of population density:
#set.seed(1003)
#
###	Check our memory footprint just to make sure:
#memory.size()
#
#
###	yaImpute Approach:
#
#start_time = proc.time()[3]
#popfit = yai(x = x_data, y = y_data, method = "randomForest", rfMode="regression", k=2)
#print(paste("Elapsed Fitting Time:", proc.time()[3] - start_time, "seconds"))
#popfit
#
#
#plot(x=impute(popfit)$y.o, y=impute(popfit)$y)
#abline(a=0, b=1, lty=2)
#
#
#yaiVarImp(popfit)
#rmsd.yai(popfit)
#yaiRFsummary(popfit)
##########



###########
###	Single core processing:
#
###	Create density estimates on 100m grid from land cover:
#
#setwd(paste(path, "Modeling Work/output", country_prefix, sep="/"))
#
###	Test the above function on a single core:
#start_time = proc.time()[3]
#
#
###	Create a raster stack of our covariates:
#covariate_stack <- stack(
#	c(
#	land_cover_v011,
#	land_cover_v040,
#	land_cover_v130,
#	land_cover_v140,
#	land_cover_v160,
#	land_cover_v190,
#	land_cover_v200,
#	land_cover_v210,
#	land_cover_v240,
#	pop_points_count,
#	pop_points_dist,
#	roads_prim_dist,
#	water_areas_perm_dist,
#	water_lines_perm_dist,
#	dem_elev,
#	dem_slope,
#	climate_zones,
#	land_cover_v011_dist,
#	land_cover_v190_dist,
#	land_cover_v240_dist,
#	land_cover_vSum,
#	census_data_zones
#	)
#)
#
#
###	Create a raster brick object based on our census data zones file to store our predictions:
#prediction_raster <- raster(census_data_zones)
###	Initialize our RasterLayer object with NAs (note that this assumes you have enough memory available to store values in memory for the entier data layer...):
#prediction_raster <- writeStart(prediction_raster, filename="predict_density.img", format="HFA", datatype="FLT4S", overwrite=TRUE)
#
#
###	NOTE: This could be accomplished by using calc() from the raster package, however it would mean predicting using our randomForest object on a per-pixel basis.. and this is slow.  So instead, I'm going to do it row by row, writing row by row...
#blocks <- blockSize(prediction_raster, chunksize=200000)
#
#pb <- tkProgressBar(title = "Predicting Population Density:", min = 0, max = blocks$n, width = 300)
#for (i in 1:blocks$n) {
#	row_data <- data.frame( getValues(covariate_stack, row=blocks$row[i], nrows=blocks$nrows[i]) )
#		
#	##	Convert field names to something more manageable and that matches our x_data:
###	Full covariate stack:
#names(row_data) <- c("PLC_011", "PLC_040", "PLC_130", "PLC_140", "PLC_160", "PLC_190", "PLC_200", "PLC_210", "PLC_240", "pp_sum", "pp_dist", "rdp_dist", "wtap_dist", "wtlp_dist", "dem_elev", "dem_slope", "DLC_011", "DLC_190", "DLC_240", "clim_zone", "PLC_Sum", "census_data_zones")
#
#
#	##	Fix factor columns:
#	row_data$clim_zone <- factor(row_data$clim_zone, levels=1:7, labels=as.character(1:7))
#
#
###	Detect if we have any NA or Inf values:
#na_present <- apply(is.na(row_data), 1, any)
#inf_present <- apply(row_data == -Inf | row_data == Inf, 1, any)
#census_data <- (is.na(row_data$census_data_zones))
#roi_subset <- (!na_present & !inf_present & !census_data)
#	
###	Create a set of predictions based on our covariates:
#predictions <- numeric(length=length(row_data[,1]))
#predictions[] <- NA
#
###	If we have data where NAs or Inf values are not present then we predict for those cells:
#if (sum(roi_subset) > 0) {
#
#	##	Convert our sums of land cover pixels to proportions:
#	##	NOTE:	Since we are using the _agg.img instead of the _focal.img
#	##		our proportion land covers are on the 100 m pixels, aggregated
#	##		from the 25 m pixels... So we derive proportion land cover by
#	##		dividing by 16 instead of the sum within the focal region...
#	#row_data[,1:9] <- row_data[,1:9] / row_data$PLC_Sum
#	row_data[,1:9] <- row_data[,1:9] / 16
#
#	##	Run predictions on data that aren't missing any data, and note that we remove the vSum column from the row_data and reorder the columns to match our x_data:
#
#	###	The retrieved data minus the vSum and census_data_zones columns is fed int
#	###		(see above):
#	#predictions[roi_subset] <- predict(popfit, newdata=row_data[roi_subset,1:(length(row_data)-2)])
#}
#
#	
#	##	Now store our predictions in our prediction
#	prediction_raster <- writeValues(prediction_raster, predictions, blocks$row[i])
#	
#	setTkProgressBar(pb, i, label=paste( round(i/blocks$n*100, 0),"% done"))
#}
#
#prediction_raster <- writeStop(prediction_raster)
#
#close(pb)
#print(paste("Elapsed Prediction Time:", proc.time()[3] - start_time, "seconds"))
#
#
###### Diagnostics:
#test_dbf <- foreign::read.dbf("predict_density_sum_by_areas.dbf")
#test_dbf$OBJECTID <- test_dbf$FID_ + 1
#merged_dbf <- merge(x=census_data@data, y=test_dbf, by.x="OBJECTID", by.y="OBJECTID", all.x=TRUE, sort=FALSE)
#merged_dbf$EST_POPD <- merged_dbf$SUM/(merged_dbf$AREA/10000)
######



#####
##	Let's parallelize the process using the snow package:
##		NOTE:  The 00.1 script used the new parallel package
##			but it seems like this is a simpler way given that we have
##			to write out our results one block at a time...

##		NOTE also:  If you've changed the predictor set then you 
##			need to change the column renaming in the cluster_predict()
##			function and the subset of the RasterLayer objects in the
##			raster brick that gets created before the process is started...

##	Create the cluster process within a function:
cluster_predict <- function(prediction_raster, ...) {
	##	Start the timer:
	start_time = proc.time()[3]

	##	Pull the cluster:
	cl <- getCluster()
	on.exit( returnCluster() )
	
	##	Determine the number of cores we're working with:
	nodes <- length(cl)

	##	Generate a set of blocks on which to process the raster
	blocks <- blockSize(prediction_raster, chunksize=200000, minblocks=nodes*4)

	pb <- tkProgressBar(title = "Predicting Population Density:", min = 0, max = blocks$n, width = 300)

	##	Pass off required libraries and data to the cluster workers:
	clusterEvalQ(cl, {
		require(raster)
		require(randomForest)
	})

	clusterExport(cl, c("popfit", "covariate_stack"))
	##	Since "blocks" only exists inside this function's environment, we need
	##		to call environment() to get our current environment rather than the 
	##		global environment assumed by clusterExport():
	clusterExport(cl, "blocks", envir=environment())

	##	Define the function that will be run on each cluster to do the predictions:
	clFun <- function (i) {

		row_data <- data.frame( getValues(covariate_stack, row=blocks$row[i], nrows=blocks$nrows[i]) )
			
		##	Convert field names to something more manageable and that matches our x_data:

		###	Full covariate stack:
		#names(row_data) <- c("PLC_011", "PLC_040", "PLC_130", "PLC_140", "PLC_160", "PLC_190", "PLC_200", "PLC_210", "PLC_240", "pp_sum", "pp_dist", "rdp_dist", "wtap_dist", "wtlp_dist", "dem_elev", "dem_slope", "DLC_011", "DLC_190", "DLC_240", "clim_zone", "PLC_Sum", "census_data_zones")
		###	Reduced covariate stack:
		#names(row_data) <- c("PLC_011", "PLC_040", "PLC_130", "PLC_140", "PLC_160", "PLC_190", "PLC_200", "PLC_210", "PLC_240", "rdp_dist", "wtap_dist", "wtlp_dist", "clim_zone", "DLC_011", "DLC_190", "DLC_240", "census_data_zones")
		###	Reduced covariate stack, Built:
		#names(row_data) <- c("PLC_011", "PLC_040", "PLC_130", "PLC_140", "PLC_160", "PLC_200", "PLC_210", "rdp_dist", "wtap_dist", "wtlp_dist", "clim_zone", "DLC_011", "PLC_Built", "DLC_Built", "census_data_zones")
		##	Reduced covariate stack, Built, Proportion "Urban" 1k:
		names(row_data) <- c("PLC_011", "PLC_040", "PLC_130", "PLC_140", "PLC_160", "PLC_200", "PLC_210", "rdp_dist", "wtap_dist", "wtlp_dist", "clim_zone", "DLC_011", "PLC_Built", "DLC_Built", "PLC1k_190", "census_data_zones")

		##	Fix factor columns:
		row_data$clim_zone <- factor(row_data$clim_zone, levels=1:7, labels=as.character(1:7))

		##	Detect if we have any NA or Inf values:
		na_present <- apply(is.na(row_data), 1, any)
		inf_present <- apply(row_data == -Inf | row_data == Inf, 1, any)
		census_data <- (is.na(row_data$census_data_zones))
		roi_subset <- (!na_present & !inf_present & !census_data)
			
		##	Create a set of predictions based on our covariates:
		predictions <- numeric(length=length(row_data[,1]))
		predictions[] <- NA

		##	If we have data where NAs or Inf values are not present then we predict for those cells:
		if (sum(roi_subset) > 0) {

			##	Convert our sums of land cover pixels to proportions:
			##	NOTE:	Since we are using the _agg.img instead of the _focal.img
			##		our proportion land covers are on the 100 m pixels, aggregated
			##		from the 25 m pixels... So we derive proportion land cover by
			##		dividing by 16 instead of the sum within the focal region...
			#row_data[,1:9] <- row_data[,1:9] / row_data$PLC_Sum
			row_data[,c(1:7,13)] <- row_data[,c(1:7,13)] / 16
			row_data[,15] <- row_data[,15] / 4048

			##	Run predictions on data that aren't missing any data, and note that we remove the vSum and/or census_data_zones column from the row_data and reorder the columns to match our x_data:

			###	The retrieved data minus the vSum and census_data_zones columns is fed int
			###		(see above):
			#predictions[roi_subset] <- predict(popfit, newdata=row_data[roi_subset,1:(length(row_data)-2)])

			##	The retrieved data minus the census_data_zones column is fed int
			##		(see above):
			predictions[roi_subset] <- predict(popfit, newdata=row_data[roi_subset,1:(length(row_data)-1)])
		}

		return(predictions)
	}

	##	Start all nodes on a prediction:
	for (i in 1:nodes) {
		sendCall(cl[[i]], clFun, i, tag=i)
	}

	##	Start the raster writer object so we can store our results as they
	##		come back from our cluster:
	setwd(paste(path, "Modeling Work/output", country_prefix, sep="/"))
	
	prediction_raster <- writeStart(prediction_raster, filename="predict_density.img", format="HFA", datatype="FLT4S", overwrite=TRUE)


	##	Create our primary cluster processing loop, recalling that we already
	##		have clusters running:
	for (i in 1:blocks$n) {
		##	Receive results from a node:
		predictions <- recvOneData(cl)

		##	Check if there was an error:
		if (!predictions$value$success) {
			stop("ERROR: Cluster barfed...")
		}

		##	Which block are we processing:
		block <- predictions$value$tag
		cat("Received block: ", block, "\n")
		flush.console()

		##	Now store our predictions in our prediction
		prediction_raster <- writeValues(prediction_raster, predictions$value$value, blocks$row[block])

		##	Check to see if we are at the end of our block list:
		ni <- nodes + i
		if (ni <= blocks$n) {
			##	And if not, send it to the cluster node that just gave us
			##		our last result...
			sendCall(cl[[predictions$node]], clFun, ni, tag=ni)
		}
	
		setTkProgressBar(pb, i, label=paste( round(i/blocks$n*100, 0),"% done"))
	}

	prediction_raster <- writeStop(prediction_raster)

	close(pb)
	print(paste("Elapsed Prediction Time:", proc.time()[3] - start_time, "seconds"))
	flush.console()

	return(prediction_raster)
}


##	Create a raster stack of our covariates and the census_data_zones file which will allow us to restrict processing to just the areas within the boundaries of our census data area (NOTE: This should be changed here to match the covariates used in the estimation of the model, as well as the renaming applied in the cluster predict function.  This will speed processing up slightly especially if used on a subset of the predictors):

##	Full covariate stack:
#covariate_stack <- stack(
#	c(
#	land_cover_v011,
#	land_cover_v040,
#	land_cover_v130,
#	land_cover_v140,
#	land_cover_v160,
#	land_cover_v190,
#	land_cover_v200,
#	land_cover_v210,
#	land_cover_v240,
#	pop_points_count,
#	pop_points_dist,
#	roads_prim_dist,
#	water_areas_perm_dist,
#	water_lines_perm_dist,
#	dem_elev,
#	dem_slope,
#	climate_zones,
#	land_cover_v011_dist,
#	land_cover_v190_dist,
#	land_cover_v240_dist,
#	land_cover_vSum,
#	census_data_zones 
#	)
#)

###	Reduced covariate stack:
#covariate_stack <- stack(
#	c(
#	land_cover_v011,
#	land_cover_v040,
#	land_cover_v130,
#	land_cover_v140,
#	land_cover_v160,
#	land_cover_v190,
#	land_cover_v200,
#	land_cover_v210,
#	land_cover_v240,
#	roads_prim_dist,
#	water_areas_perm_dist,
#	water_lines_perm_dist,
#	climate_zones,
#	land_cover_v011_dist,
#	land_cover_v190_dist,
#	land_cover_v240_dist,
#	census_data_zones
#	)
#)

###	Reduced covariate stack, Built:
#covariate_stack <- stack(
#	c(
#	land_cover_v011,
#	land_cover_v040,
#	land_cover_v130,
#	land_cover_v140,
#	land_cover_v160,
#	land_cover_v200,
#	land_cover_v210,
#	roads_prim_dist,
#	water_areas_perm_dist,
#	water_lines_perm_dist,
#	climate_zones,
#	land_cover_v011_dist,
#	land_cover_v190v240,
#	land_cover_v190v240_dist,
#	census_data_zones
#	)
#)

##	Reduced covariate stack, Built:
covariate_stack <- stack(
	c(
	land_cover_v011,
	land_cover_v040,
	land_cover_v130,
	land_cover_v140,
	land_cover_v160,
	land_cover_v200,
	land_cover_v210,
	roads_prim_dist,
	water_areas_perm_dist,
	water_lines_perm_dist,
	climate_zones,
	land_cover_v011_dist,
	land_cover_v190v240,
	land_cover_v190v240_dist,
	land_cover_v190_f,
	census_data_zones
	)
)

##	Create a raster object based on our census data zones file to store our predictions:
prediction_raster <- raster(census_data_zones)


##	Set the extent of processing (you can crop here if you want):
my_extent <- extent(census_data_zones)
#my_extent <- extent(census_data_zones, 7000, 8000, 7000, 8000)
#my_extent <- extent(census_data_zones, 6000, 9000, 6000, 9000)
#plot(census_data_zones, ext=my_extent)
prediction_raster <- crop(prediction_raster, my_extent)
covariate_stack <- crop(covariate_stack, my_extent)

beginCluster()
prediction_raster <- cluster_predict(prediction_raster)
endCluster()

#Received block:  833 
#[1] "Elapsed Prediction Time: 16625.25 seconds"

##	4.61 hours for a raster of 12492 x 13236 or 165344112 cells on a 
##		dual core laptop.  It	takes just under twice as long to run on a 
##		single core, with perhaps	a half hour of overhead tacked on.  More
##		complex randomForest models will take longer to run (up to 6.5 hours
##		dual core on the LUECI desktops).


##	Save the workspace:
setwd(paste(path, "Modeling Work/output", country_prefix, sep="/"))

#save.image(file="00.2_image.RData")
load(file="00.2_image.RData")
#####
