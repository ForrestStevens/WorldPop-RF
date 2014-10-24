#####
##	BEGIN:	Load configuration options

##	This section can be safely commented as long as you have already opened
##		and run the Metadata.R file that should be located within the 
##		RF/data/*/ folder where * represents the country name for which the
##		RF folder is located.

##	NOTE:  It's critical that you realize that anything you specify here 
##		will be overridden by the values read in from the Metadata.r file
##		so the paths may be incorrect if specified incorrectly inside that
##		file!

##	Configure the country abbreviation and name:
country <- "KHM_002"

##	This should be set to the folder *containing* the "RF" folder structure:
root_path <- "D:/Documents/Graduate School/Research/Population/Data/"
project_path <- paste(root_path, "RF/data/", country, "/", sep="")

##	Load the metadata from the file created in the country's /data folder:
source(paste(project_path, "Metadata.r", sep=""))

##	END:	Load configuration options
#####



#####
##	BEGIN:	Validation configuration options

##	NOTE:	This script assumes that you have prepared output at a coarse
##		level with which to compare to census data in a country as specified
##		here.
country_comp <- "KHM"

##	END:	Validation configuration options
#####




#####
##	NOTICE:  In practice nothing below this line should need to be regularly
##		edited unless very specific details need to be changed about the 
##		modeling process (e.g. speciyfing an extent or something detailed
##		about the RandomForest modeling).
#####



#####
##	BEGIN:	Package loading and configuration

require(rgdal)
require(raster)
require(RPyGeo)
require(foreign)
require(ggplot2)


##	Parameters and defaults:


##	Setup output paths for the coarse level and project paths for the
##		finest level census data:

##	Default coarse level output (for automated output from my scripts):
output_file <- paste(output_path, country, "_popmap_v", version, ".tif", sep="")


##	Default census file for the finest level available (folder is the country
##		code listed without a level of aggregation):

##	This should be set to the folder *containing* the "RF" folder structure:
project_path_comp <- paste(root_path, "RF/data/", country_comp, "/", sep="")


dataset_name <- list.files(paste(project_path_comp, "Census/", sep=""),  "shp$")
census_file <- paste(project_path_comp, "Census/", dataset_name, sep="")


##	END:	Package loading and configuration
#####


#####
##	BEGIN:	Internal function declaration


##	END:	Internal function declaration
#####



#####
##	BEGIN:	Metadata setup and compilation


##	Clean up temporary variables:


##	END:	Metadata setup and compilation
#####



#####
##	BEGIN:	Data loading and data structure setup

##	If our covariate data already exist load it, if not create it:
setwd(paste(project_path, "Census", sep=""))


workspace_dir <- "C:/tmp"

zonal_raster <- paste(workspace_dir, "/census_zones.tif", sep="")
zonal_points <- paste(workspace_dir, "/census_zones_points.shp", sep="")
zonal_raster_stats <- paste(workspace_dir, "/zonal.tif", sep="")
zonal_points_extract <- paste(workspace_dir, "/census_zones_points_extract.shp", sep="")


##	We don't assume that we've already processed the census file and
##		have the zonal pieces created.  So we'll create them on the fly
##		here:

##	NOTE:	It appears I've found the culprit with the ArcGIS/Geoprocessing
##		Zonal Statistics and other table writing functions to inconsistently
##		crash with a read/write locking error.  If you turn off Microsoft
##		Security Essentials' "Real Time Protection" then the writes complete
##		just fine.  So I've switched back to using the custom zonal statistics
##		as the process is *much* faster than the zonal() approach:

rpygeo_env <- rpygeo.build.env(extensions="Spatial", python.path="C:/Python26/ArcGIS10.0", python.command="python.exe", workspace=workspace_dir, overwriteoutput=1)

rpygeo.geoprocessor(
	paste(
		##	This generates a zonal raster and then uses a converted feature to points (inside) file to extract the zonal data (note that with multiple lines you need to add the four space indent):
		"import arcpy\n",
		"    arcpy.CheckOutExtension(\"Spatial\")\n",
		"    from arcpy.sa import *\n",
		"    arcpy.FeatureToRaster_conversion(\"", census_file, "\", \"ADMINID\", \"", zonal_raster, "\", 0.0008333)\n",
		"    gp.ZonalStatistics_sa(\"", zonal_raster, "\", \"Value\", \"", output_file, "\" ,\"", zonal_raster_stats, "\", \"SUM\", \"DATA\")\n",
		"    arcpy.FeatureToPoint_management(\"", census_file, "\", \"", zonal_points, "\", \"INSIDE\")\n",
		"    gp.ExtractValuesToPoints_sa(\"", zonal_points, "\", \"", zonal_raster_stats, "\", \"", zonal_points_extract, "\", \"NONE\", \"VALUE_ONLY\")",
	sep=""),
	env=rpygeo_env, 
	add.gp=FALSE, 
	clean.up=FALSE, 
	working.directory=workspace_dir
)

##	Check to see if there's an error and bail:
if (file.exists(paste(workspace_dir, "/rpygeo.msg", sep=""))) {
	stop("ERROR:  Error in geoprocessing the zonal statistics... Check the .msg file in the geoprocessing working directory!")
	flush.console()
} else {
	##	Now load up the DBF and return it as a data.frame:
	output_dbf <- foreign::read.dbf(sub(".shp", ".dbf", zonal_points_extract))
}


##	END:	Data loading and data structure setup
#####



#####
##	BEGIN:	Validation of data and output:
observed <- output_dbf$ADMINPOP
predicted <- output_dbf$RASTERVALU

##	Calculate RMSE:
rmse <- sqrt(sum((observed - predicted)^2) / length(observed))
rmse

##	Calculate %RMSE (percentage of the mean census block size):
pct_rmse <- rmse/mean(observed)*100
pct_rmse

##	Calculate MAE:
mae <- sum(abs(observed - predicted)) / length(observed)
mae

##	Scatterplot of data:
#qplot(predicted, observed)

png(file=paste(output_path, "predicted_vs_observed_v", pop_version, ".png", sep=""))
plot(y=predicted, x=observed, col=rgb(0,0,0,0.2), xlim=c(min(c(observed, predicted)), max(c(observed, predicted))), ylim=c(min(c(observed, predicted)), max(c(observed, predicted))), xlab="Observed", ylab="Predicted", pch=16, cex=0.7, main=paste(country, "_", aggregation_level, sub_data, sep=""))
abline(a=0, b=1, lty=2, col="darkgrey")

text(y=max(c(observed, predicted))/10*4, x=max(c(observed, predicted))/10*5.5, paste("RMSE: ", format(rmse, nsmall=2), sep=""), pos=4)
text(y=max(c(observed, predicted))/10*3, x=max(c(observed, predicted))/10*5.5, paste("%RMSE: ", format(pct_rmse, nsmall=2), sep=""), pos=4)
text(y=max(c(observed, predicted))/10*2, x=max(c(observed, predicted))/10*5.5, paste("MAE: ", format(mae, nsmall=2), sep=""), pos=4)
dev.off()


##	Repeat the plot for visualizing in R:
plot(y=predicted, x=observed, col=rgb(0,0,0,0.2), xlim=c(min(c(observed, predicted)), max(c(observed, predicted))), ylim=c(min(c(observed, predicted)), max(c(observed, predicted))), xlab="Observed", ylab="Predicted", pch=16, cex=0.7, main=paste(country, "_", aggregation_level, sub_data, sep=""))
abline(a=0, b=1, lty=2, col="darkgrey")

text(y=max(c(observed, predicted))/10*4, x=max(c(observed, predicted))/10*5.5, paste("RMSE: ", format(rmse, nsmall=2), sep=""), pos=4)
text(y=max(c(observed, predicted))/10*3, x=max(c(observed, predicted))/10*5.5, paste("%RMSE: ", format(pct_rmse, nsmall=2), sep=""), pos=4)
text(y=max(c(observed, predicted))/10*2, x=max(c(observed, predicted))/10*5.5, paste("MAE: ", format(mae, nsmall=2), sep=""), pos=4)


##	END:	Validation of data and output:
#####
