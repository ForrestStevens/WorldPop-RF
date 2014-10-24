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


##	Parse main configuration file, which will set the country and root_path
##		variables:
source("01.0 - Configuration.py.r")


##	Load the metadata from the file created in the country's /data folder:
project_path <- paste(root_path, "/data/", country, "/", sep="")
source(paste(project_path, "Metadata.r", sep=""))

##	END:	Load configuration options
#####



#####
##	BEGIN:	RandomForest configuration

##	Configuration options for RandomForest modeling and prediction:

##	NOTE:  The following were moved to the Metadata.r file for per-country
##		configuration and reporting purposes:

##  If we are using a set of covariates from another country set the 
##		fixed_set variable to specify which will then be used to fix the 
##		covariate list to an existing randomForest object, otherwise, 
##		use the full set from the covariate metadata if it is NULL.
##
##    Note that the fixed_set flag also changes the way that the 
##		randomForest object is created below by eliminating the variable 
##		elimination routine:
#fixed_set <- "VNM"
#fixed_set <- c("VNM", "KHM")
#fixed_set <- NULL
if (!exists("fixed_set")) {
	fixed_set <- NULL
}


##	END:	RandomForest configuration
#####



#####
##	NOTICE:  In practice nothing below this line should need to be regularly
##		edited unless very specific details need to be changed about the 
##		modeling process (e.g. speciyfing an extent or something detailed
##		about the RandomForest modeling).
#####



#####
##	BEGIN: Package loading and fixed parameter configuration

require(rgdal)
require(raster)
require(RPyGeo)
require(randomForest)
require(quantregForest)
require(foreign)
require(snow)
##	NOTE: The tcltk library is used to create a progress bar, but the
##		raster package has a progress bar mechanism built in in recent versions.
##		Look up the pbCreate(), pbStep(), and pbClose() functions in the raster
##		package...
require(tcltk)


##	Fixed parameters and defaults:
proj4str_gcs_wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


##	Set the path to the Python location which has access to ArcGIS 
##		Geoprocessing facilities.  As long as you're running RStudio using
##		the batch file to start it the ARCPY environment variable should
##		contain the appropriate path:
#python_path <- "C:/Python26/ArcGIS10.0"
python_path <- Sys.getenv("ARCPY")

if (python_path == "") {
	stop("ERROR:  There's an error somewhere in your configuration!  The system path to Python does not exist.  The RF software expects to have the ARCPY= environment setting set either by the system or the RStudio batch file to the location of ArcGIS' version of Python, by default, C:/Python27/ArcGIS10.1 for example for ArcGIS 10.1.  Double check that ARCPY is set in the environment settings and/or in the RStudio batch file!")
}

##	END: Package loading and fixed parameter configuration
#####



#####
##  BEGIN:	Internal function declaration

##	Transform Y:
transY <- function(x, inverse=FALSE) {
  ###	The default is to not transform Y:
  #return(x)
  
  ##	The next one we'll try is log transforming:
  if (!inverse) {
    return( log(x) )
  } else {
    return( exp(x) )
  }
}

##	NOTE:	This function is no longer needed because of the inconsistency 
##		with which the geoprocessing zonal statistics tool runs on large
##		datasets.  I replaced it using the slower, but more foolproof raster
##		package zonal() function in the code:

##	NOTE:	It appears I've found the culprit with the ArcGIS/Geoprocessing
##		Zonal Statistics and other table writing functions to inconsistently
##		crash with a read/write locking error.  If you turn off Microsoft
##		Security Essentials' "Real Time Protection" then the writes complete
##		just fine.  So I've switched back to using the custom zonal statistics
##		as the process is *much* faster than the zonal() approach:

arcgis_zonal <- function(raster_file, census_file, mask_file="", stat="MEAN") {
	workspace_dir <- "C:/tmp"

	rpygeo_env <- rpygeo.build.env(extensions="Spatial", python.path=python_path, python.command="python.exe", workspace=workspace_dir, overwriteoutput=1)

	##	Basic fixes for accomodating the R zonal() arguments:
	stat = toupper(stat)
	if (stat=="MODAL") stat = "MAJORITY"

	##	If we are using a mask file we need to first multiply our mask through our output file:
	if (mask_file != "") {
		rpygeo.geoprocessor(
			paste(
				"import arcpy\n",
				"    arcpy.CheckOutExtension(\"Spatial\")\n",
				"    from arcpy.sa import *\n",
				"\n",
				"    mask_raster = Raster(\"", mask_file, "\")\n",
				"    data_raster = Raster(\"", raster_file, "\")\n",
				"    out_raster = SetNull(mask_raster == 1, data_raster)\n",
				"    out_raster.save(\"", workspace_dir, "/out_raster.tif\")\n",
				"    out_raster = None",
			sep=""),
			env=rpygeo_env, 
			add.gp=FALSE, 
			clean.up=FALSE, 
			working.directory=workspace_dir
		)
		
		##	Set our raster_file to our masked file:
		raster_file <- paste(workspace_dir, "out_raster.tif", sep="/")
		
		##	Check to see if there's an error and bail:
		if (file.exists(paste(workspace_dir, "/rpygeo.msg", sep=""))) {
			stop("ERROR:  Error in geoprocessing the zonal statistics... Check the .msg file in the geoprocessing working directory!")
			flush.console()

			return(NULL)
		}
	}
	
	
	rpygeo.geoprocessor(
		paste(
			##	Use the first if you're using a shapefile with an ADMINID field, the second if you've already converted to a zonal raster:
			#"gp.ZonalStatisticsAsTable_sa(\"", census_file, "\",\"ADMINID\",\"", raster_file,"\",\"", workspace_dir, "/my_zonalstats_output.dbf\",\"DATA\",\"", stat, "\")", 

			##	This option still works inconsistently in ArcGIS 10.0, so we've transitioned to using the third option below:
			#"gp.ZonalStatisticsAsTable_sa(\"", census_file, "\",\"Value\",\"", raster_file,"\",\"", workspace_dir, "/my_zonalstats_output.dbf\",\"DATA\",\"", stat, "\")", 

			##	This generates a zonal raster and then uses a converted feature to points (inside) file to extract the zonal data (note that with multiple lines you need to add the four space indent):
			"gp.ZonalStatistics_sa(\"", census_file, "\", \"Value\",\"", raster_file, "\", \"", workspace_dir, "/zonal.tif\", \"", stat, "\", \"DATA\")\n",
			"    gp.ExtractValuesToPoints_sa(\"", sub("\\_zones..*", "_points.shp", census_file), "\", \"", workspace_dir, "/zonal.tif\", \"", workspace_dir, "/zonal_points.shp\", \"NONE\", \"VALUE_ONLY\")",
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

		return(NULL)
	}

	##	Now load up the DBF and return it as a data.frame:
	#output_dbf <- foreign::read.dbf(paste(workspace_dir, "/my_zonalstats_output.dbf", sep=""))
	output_dbf <- foreign::read.dbf(paste(workspace_dir, "/zonal_points.dbf", sep=""))

	##	Formatted output to match the R raster zonal() function:
	##		Used for the first two zonal options above:
	#return(output_dbf[,c(1,4)])
	##		Used for the points extraction method:
	return(output_dbf[,c("ADMINID", "RASTERVALU")])
}

##	END:	Internal function declaration
#####



#####
##	BEGIN:	Data loading and data structure setup

##	Load the compiled covariate data as constructed from the Metadata.r file
##		in the "Data Preparation, R.r" script:
load(file=paste(output_path_tmp, "covariates.RData", sep=""))


##	If our covariate data already exist load it, if not create it:
setwd(paste(project_path, "Census/Derived", sep=""))

##	Read in our zonal raster based on ADMINID values, this is used
##		for zonal statistics and as a template for predictions:
zonal_raster_path <- paste(project_path, "Census/Derived/census_zones.tif", sep="")
zonal_raster <- raster(zonal_raster_path)


start_time = proc.time()[3]

if (file.exists("census_covariates.shp")) {
	census_data <- readOGR(".", "census_covariates")
} else {
  ##  You can set doit to FALSE if you want to start processing the 
  ##    covariates from somewhere in the middle (e.g. if the zonal stats
  ##    barfed).  If doit is FALSE then assume that census_data already
  ##    is partially processed and don't load it again.)  You should set the
	##		doit_var_name argument to the last known successfully processed
	##		variable:
  doit <- TRUE
  #doit <- FALSE
	#doit_var_name <- "roa_dst"


  if (doit) {
	  ##	Read in our census data with final population counts and area field, F_AREA (in m^2):
	  census_data <- readOGR(".", "census_area")
  }

	##	Recalculate our areas into hectares and add them to the census data:
	census_data$AREA_HA <- census_data$F_AREA / 10000

	##	Finally calculate our population density in people per hectare for use as a
	##		covariate variable:
	census_data$POPD_PPHA <- census_data$ADMINPOP / census_data$AREA_HA

	##	Set our mask path based on water areas so that water areas don't
	##		count in our summarizations by admin unit.  If you want this
	##		turned off just set it to "":
	#mask_path <- ""
	mask_path <- covariates[["lan_cls210"]]$path
	
	##	For each dataset in our covariate data structure, we need to run our custom Geoprocessing zonal statistics over the datasets:
	for ( i in 1:length(covariates) ) {
		var_name <- names(covariates[i])
		
		covariate <- covariates[[i]]
		dataset_summary <- covariate$dataset_summary
		
		raster_path <- covariate$path
		dataset_raster <- raster(raster_path)

		##	If this is a 'modal' summary and there exists a class
		##		summary for the same value already calculated in the output
		##		data.frame, then just calculate whether the majority will be 1-0 
		##		based on the proportion rather than incurring the overhead of 
		##		the zonal statistics.  Note that for the zonal() the function
		##		is modal(), but if using arcgis_zonal() then the function
		##		is MAJORITY:
		
		##	If you need to subset or restart you can set doit=FALSE above
		##		and then set the doit_var to a particular var_name to start at 
		##		particular place:
		if (!doit) {
			if (var_name == doit_var_name) {
				doit = TRUE
			}
		}
		if (doit) {
			if (dataset_summary == "modal" & gsub("cls", "prp", var_name) %in% names(census_data)) {
				census_data@data[[var_name]] <- as.numeric(census_data@data[[ gsub("cls", "prp", var_name) ]] > 0.5)
			} else {
				#output_stats <- zonal(dataset_raster, zonal_raster, stat=dataset_summary)
				output_stats <- arcgis_zonal(raster_path, zonal_raster_path, mask_file=mask_path, stat=dataset_summary)
				##	Store the output summary column in our census_data:
				##	NOTE:	We can no longer assume that the ADMINID field is sorted, 
				##		so we have to do a merge and use the resulting output to match:
				names(output_stats) <- c("ADMINID", dataset_summary)
				output_stats <- merge(as.data.frame(census_data)['ADMINID'], output_stats, by="ADMINID", sort=FALSE)
				
				census_data@data[[var_name]] <- output_stats[,2]
				
			}
		}

		print(paste("Finished summarizing", var_name, "by zone..."))
		flush.console()
	}

	##	Cleanup:
	#rm(list=c("i", "dataset_raster", "output_stats", "var_name", "covariate", "dataset_summary", "raster_path"))


	##	Write out our shapefile with covariates:
	writeOGR(census_data, ".", paste("census_covariates", sep=""), driver="ESRI Shapefile", check_exists=TRUE, overwrite_layer=TRUE)
}

print(paste("Elapsed Processing Time:", proc.time()[3] - start_time, "seconds"))

##	END:	Data loading and data structure setup
#####



#####
##	BEGIN:	Estimate RandomForest model
setwd(output_path)

#####
##	Set up response and covariate dataframes for the random forest modeling:
y_data <- census_data@data$POPD_PPHA


if (!is.null(fixed_set)) {
  country_old <- fixed_set[1]
  output_path_tmp <- paste(root_path, "/output/", country_old, "/tmp/", sep="")
  
	##	Load the first listed, existing country and rename its popfit_final:
	load(file=paste(output_path_tmp, "popfit_final.RData", sep=""))
  popfit_final_old <- popfit_final

  fixed_predictors <- row.names(importance(popfit_final_old))
  
  output_path_tmp <- paste(root_path, "/output/", country, "/tmp/", sep="")
} else {
  fixed_predictors <- names(covariates)
}


##  Full covariate set:
x_data <- census_data@data[fixed_predictors]


##	Convert to factors where appropriate:
for (var_name in names(x_data)) {
	if (grepl("cls", var_name, perl=T)) {
		x_data[[var_name]] <- factor(x_data[[var_name]], levels=0:10, labels=as.character(0:10))
	}
}


##	NOTE: Rarely this needs to be applied because of significatly strange
##		relationships between proportion rural and other types of interactions
##		resulting in odd partial plot relationships:
##	Subset data to exclude all proportion covariates:

###	Remove all proportion data:
x_data <- x_data[,!grepl("prp", names(x_data))]
##	Remove only proportion rural data:
#x_data <- x_data[,!(names(x_data) == "lan_prp240")]

###	Remove all class data:
#x_data <- x_data[,!grepl("cls", names(x_data))]

###	Remove rural/urban built data, leaving only BLT:
x_data <- x_data[,!grepl("190", names(x_data))]
x_data <- x_data[,!grepl("240", names(x_data))]

###  Remove health facility data:
#x_data <- x_data[,!grepl("^hea_", names(x_data))]

###  Remove AfriAsiaPop gridx:
#x_data <- x_data[,!(names(x_data) == "afr")]

###  Remove elevation:
#x_data <- x_data[,!(names(x_data) == "ele")]

###  Remove NPP:
#x_data <- x_data[,!(names(x_data) == "npp")]


###	Remove all populated place data:
#x_data <- x_data[,!grepl("pop", names(x_data))]

###	Subset x_data columns to land cover alone:
#x_data <- x_data[,c(1:28)]
#x_data <- x_data[,grepl("lan_", names(x_data), perl=T)]

###	Subset x_data to remove NAs:
indexX <- complete.cases(x_data)


##	Subset y_data to remove any of interest:
#indexY <- y_data < 10 		##	Subset data to remove outliers...
#indexY <- y_data < 100 	##	Subset data to remove outliers...
indexY <- y_data > 0 			##	Subset to remove zeroes, good for 
													##		log transformation, since it results 
													##		in -Inf values...
#indexY <- y_data > -1 		##	Don't subset data...
#indexY <- y_data < transY(4, inverse=TRUE)


##	Subset data according to indices:
y_data <- y_data[indexX & indexY]
x_data <- x_data[indexX & indexY,]


##	Create a sampling factor that allows us to stratify our sample
##		during randomForest estimation, that captures as many large 
##		values as possible:
##
#y_sample <- as.factor(round((y_data/max(y_data)*10)))
y_sample <- as.factor(1)


##	Transform y_data (as defined in the transY() function:
y_data <- transY(y_data)



##	Random Forest regression of population density:
set.seed(2002)

##	Check our memory footprint just to make sure:
memory.size()


##	Now we are going tune our randomForest population density regression: 
start_time = proc.time()[3]
#init_popfit = tuneRF(x=x_data, y=y_data, plot=TRUE, mtryStart=length(x_data)/3, ntreeTry=length(y_data)/20, improve=0.0001, stepFactor=1.20, trace=TRUE, doBest=TRUE, nodesize=length(y_data)/1000, na.action=na.omit, importance=TRUE, proximity=FALSE, sampsize=length(y_data), replace=TRUE) 
init_popfit = tuneRF(x=x_data, y=y_data, plot=TRUE, mtryStart=length(x_data)/3, ntreeTry=length(y_data)/20, improve=0.0001, stepFactor=1.20, trace=TRUE, doBest=TRUE, nodesize=length(y_data)/1000, na.action=na.omit, importance=TRUE, proximity=FALSE, sampsize=min(c(length(y_data), 1000)), replace=TRUE) 
print(init_popfit)
print(paste("Elapsed Fitting Time:", proc.time()[3] - start_time, "seconds"))


##	Save off our init_popfit object for this set of data:
save(init_popfit, file=paste(output_path_tmp, "init_popfit.RData", sep=""))
#load(file=paste(output_path_tmp, "init_popfit.RData", sep=""))


if (is.null(fixed_set)) {
  ##	Now we will optimize the model by iteratively removing any 
  ##		covariates with negative increases in node purity:
  
  ##	Get list of covariates that have an importance score greater than 0:
  importance_scores <- importance(init_popfit)[order(importance(init_popfit)[,1], decreasing=TRUE),]
  pos_importance <- rownames(importance_scores)[importance_scores[,1] > 0]
  
  while (length(pos_importance) < length(importance_scores[,1])) {
  	##	Subset our x_data to just those columns having positive scores:
  	x_data <- x_data[pos_importance]
  	#popfit = tuneRF(x=x_data, y=y_data, plot=TRUE, mtryStart=length(x_data)/3, ntreeTry=length(y_data)/20, improve=0.0001, stepFactor=1.20, trace=TRUE, doBest=TRUE, nodesize=length(y_data)/1000, na.action=na.omit, importance=TRUE, proximity=FALSE, sampsize=length(y_data), replace=TRUE) 
  	popfit = tuneRF(x=x_data, y=y_data, plot=TRUE, mtryStart=length(x_data)/3, ntreeTry=length(y_data)/20, improve=0.0001, stepFactor=1.20, trace=TRUE, doBest=TRUE, nodesize=length(y_data)/1000, na.action=na.omit, importance=TRUE, proximity=FALSE, sampsize=min(c(length(y_data), 1000)), replace=TRUE) 
  	
  	##	Re-check importance scores:
  	importance_scores <- importance(popfit)[order(importance(popfit)[,1], decreasing=TRUE),]
  	pos_importance <- rownames(importance_scores)[importance_scores[,1] > 0]
  	print(popfit)
  }
  print(paste("Elapsed Fitting Time:", proc.time()[3] - start_time, "seconds"))

} else {

  popfit = randomForest(x=x_data, y=y_data, mtry=popfit_final_old$mtry, ntree=popfit_final_old$ntree, nodesize=length(y_data)/1000, importance=TRUE, proximity=TRUE)
  print(popfit)
  
}

##	Save off our popfit object for this set of data:
save(popfit, file=paste(output_path_tmp, "popfit.RData", sep=""))
#load(file=paste(output_path_tmp, "popfit.RData", sep=""))


###	Check our diagnostics:
###	Recall that running predict on the randomForest object will return OOB predictions:
#
#summary(popfit)
#plot(popfit)
#
#predict(popfit)
#importance(popfit)
#
###	Get the variable names of the top 20 predictors:
#names(importance(popfit)[,"%IncMSE"][order(importance(popfit)[,"%IncMSE"], decreasing=TRUE)])[1:20]
#
###	Sort all covariates by they $IncMSE:
#importance(popfit)[order(importance(popfit)[,1], decreasing=TRUE),]


#varImpPlot(popfit)
#varUsed(popfit)
#summary(treesize(popfit))
#plot(popfit)

# ##	Partial Plots:
# for (var_name in names(x_data)) {
# 	eval( parse( text=
# 		paste(
# 			"partialPlot(x=popfit, pred.data=x_data, x.var=\"",
# 			as.character(var_name),
# 			"\")"
# 		, sep="")
# 	))
# }


##	For continuous regression, plot observed vs. predicted:
plot(x=y_data, y=predict(popfit), ylim=c(min(y_data),max(y_data)), xlim=c(min(y_data),max(y_data)))
abline(a=0, b=1, lty=2)

###	For continuous regression, plot residuals vs. observed:
#plot(x=y_data, y=(y_data - predict(popfit)), xlim=c(0,max(y_data)))
#abline(a=0, b=0, lty=2)
#
###	For continuous regression, plot residuals vs. fitted:
#plot(x=predict(popfit), y=(y_data - predict(popfit)), xlim=c(0,max(y_data)))
#abline(a=0, b=0, lty=2)

varImpPlot(popfit)


###	Note that the predict(popfit) give predictions using OOB estimates
###		that are stored along with the randomForest object.  If we want
###		predictions from the non-OOB estimates we can run predict() using 
###		the existing x_data as the newdata= object.  This will also allow
###		us to extract the SD of the estimates from all forests, giving us
###		a pseudo-prediction uncertainty for each object/pixel we predict:
#x_data_new <- as.data.frame(census_data[names(popfit$forest$xlevels)])
#
###	Convert to factors where appropriate:
#for (var_name in names(x_data_new)) {
#	if (grepl("cls", var_name, perl=T)) {
#		x_data_new[[var_name]] <- factor(x_data_new[[var_name]], levels=0:10, labels=as.character(0:10))
#	}
#}
#
#predictions <- predict(popfit, newdata=x_data_new, predict.all=TRUE)
#
###	The wrinkle is that these predictions are the mean of the 
###		non-backtransformed predictions from each of the trees in the
###		forest.  This is problematic since the back-transform of the mean
###		is not the same as the mean of the back-transformed predictions, 
###		but this is probably what we want because taking the mean or SD
###		of the back-transformed values predicted by each tree is difficult 
###		to rationalize and based on some preliminary testing actually fits
###		the aggregated data worse than taking the mean of the individual
###		tree predictions and backtransforming the mean to use as the 
###		predicted value.
#
###	Untransformed mean of predictions for first census unit:
#predictions$aggregate[1]
###		is the same as:
#mean(predictions$individual[1,])
###		However, since we are trying to get the value from a transformed
###		response variable, these are equivalent:
#transY(predictions$aggregate[1], inverse=TRUE)
#transY(mean(predictions$individual[1,]), inverse=TRUE)
###		But not to this (which is not as accurate for approximately log-normal
###			variables, so we'll stick with the first way):
##mean(transY(predictions$individual[1,], inverse=TRUE))
#
###	To generate mean and CV for predictions across all forests for each unit:
##census_data$rf_pred <- apply(transY(predictions$individual, inverse=TRUE), MARGIN=1, mean)
###		is not equivalent to, and a more biased way than the back 
###		transformation of the mean of the predictions, which is actually 
###		the median of a log-normal backtransformed prediction distribution, 
###		which is probably what we really want given that it's more accurate:
#census_data$rf_pred <- transY( apply(predictions$individual, MARGIN=1, mean), inverse=TRUE)
#
###	Calculate residuals:
#census_data$rf_res <- census_data[["POPD_PPHA"]] - census_data[["rf_pred"]]
#
###	This is completely wrong for a non-linear, log transformation, as 
###		the back-transform of the SD of the non-back-transformed predictions 
###		is meaningless:
##census_data$rf_sd <- transY(apply(predictions$individual, MARGIN=1, sd), inverse=TRUE)
###		So instead we could back-transform and then calculate the SD:
##census_data$rf_sd <- apply(transY(predictions$individual, inverse=TRUE), MARGIN=1, sd)
##census_data$rf_cv <- (census_data[["rf_sd"]] / census_data[["rf_pred"]]) * 100
###			but the problem here is that the SD or CV for a highly skewed, log-
###			normal distribution is not a useful measure of dispersion.  So
###			instead we'll just calculate the SD of the log transformed 
###			predictions, and the CV of that value:
#census_data$rf_sd <- apply(predictions$individual, MARGIN=1, sd)
#
###		NOTE: The CV of a log normal distribution is not calculated as the
###			SD/Mean, because the distribution is not ratio-scale.  Therefore 
###			this is INCORRECT:
##census_data$rf_cv <- apply(predictions$individual, MARGIN=1, sd) / apply(predictions$individual, MARGIN=1, mean)
###			And this does not have any meaning on the interval scale so we 
###			have to use the CV on the original scale using this relationship:
###				(http://www.stata.com/statalist/archive/2003-06/msg00508.html)
#census_data$rf_cv <- sqrt(exp( (apply(predictions$individual, MARGIN=1, sd))^2 ) - 1)
###		Or alternatively, could calculate the "Geometric CV", but this seems
###			to be frowned on (ala Wikipedia's entry on CV):
###			(http://www.sportsci.org/resource/stats/logtrans.html)
##census_data$rf_cv <- exp( apply(predictions$individual, MARGIN=1, sd) ) - 1
#
###	Plot means and CV for each census unit, with no border (col=NA) around
###		polygons:
#spplot(census_data, zcol=c("POPD_PPHA"), col=NA)
#spplot(census_data, zcol=c("rf_pred"), col=NA)
#spplot(census_data, zcol=c("rf_res"), col=NA)
#spplot(census_data, zcol=c("rf_sd"), col=NA)
#spplot(census_data, zcol=c("rf_cv"), col=NA)


##	Another alternative is to use Quantile Regression Forests to generate
##		prediction intervals.  We'll fit a quantile regression using
##		the tuning parameters pulled from the popfit object above:
set.seed(2010)
popfit_final = randomForest(x=x_data, y=y_data, mtry=popfit$mtry, ntree=popfit$ntree, nodesize=length(y_data)/1000, importance=TRUE, proximity=TRUE)
set.seed(2010)
popfit_quant = quantregForest(x=x_data, y=y_data, mtry=popfit$mtry, ntree=popfit$ntree, nodesize=length(y_data)/1000)

#popfit_quant
#plot(popfit_quant)
#varImpPlot(popfit_quant)
#mean(popfit_quant$rsq)
#popfit_final
#
###	Plot predicted mean log(density) to observed (standard behavior of 
###		randomForest):
#plot(x=y_data, y=predict(popfit_final), ylim=c(min(y_data),max(y_data)), xlim=c(min(y_data),max(y_data)))
#abline(a=0, b=1, lty=2)
#
###	Plot predicted median log(density) to observed:
#plot(x=y_data, y=predict(popfit_quant)[,2], ylim=c(min(y_data),max(y_data)), xlim=c(min(y_data),max(y_data)))
#abline(a=0, b=1, lty=2)
#
###	Compares the means and medians of our predictions (remember these are log transformed):
#plot(x=predict(popfit_final), y=predict(popfit_quant)[,2], ylim=c(min(y_data),max(y_data)), xlim=c(min(y_data),max(y_data)))
#abline(a=0, b=1, lty=2)
#
#randomForest:::print.randomForest(popfit_quant)
#
###	Generate predicted values for new data and their 5% and 95% quantiles:
#predict(popfit_quant, newdata=x_data)
#
###	Predict our median and 90% lower and upper quantiles:
#predict(popfit_quant, newdata=x_data, quantiles=c(.05,0.5,.95))

##	Save off our popfit object for this set of data:
save(popfit_final, file=paste(output_path_tmp, "popfit_final.RData", sep=""))
#load(file=paste(output_path_tmp, "popfit_final.RData", sep=""))
save(popfit_quant, file=paste(output_path_tmp, "popfit_quant.RData", sep=""))
#load(file=paste(output_path_tmp, "popfit_quant.RData", sep=""))


##########
## Set the fixed_set to existing countries if you are using an existing
##    set of randomForest objects to predict from:
if (!is.null(fixed_set)) {
	if ((length(fixed_set) != 1) | !(country %in% fixed_set)) {
		##	Check to see if the country we are processing is included in our 
		##		fixed_set and if it is not then pull the first country from
		##		the set to use as our starting RF model to combine with:
		if (country %in% fixed_set) {
			index <- 1

			popfit_final_combined <- popfit_final
			popfit_quant_combined <- popfit_quant
		} else {
			index <- 2
			
			country_old <- fixed_set[1] 
			output_path_tmp <- paste(root_path, "/output/", country_old, "/tmp/", sep="")
			load(file=paste(output_path_tmp, "popfit_final.RData", sep=""))
			popfit_final_combined <- popfit_final
			load(file=paste(output_path_tmp, "popfit_quant.RData", sep=""))
			popfit_quant_combined <- popfit_quant
		}

		while (index <= length(fixed_set)) {
			## Load randomForest objects to combine:
			country_old <- fixed_set[index] 

			if (country_old != country) {
				output_path_tmp <- paste(root_path, "/output/", country_old, "/tmp/", sep="")
				load(file=paste(output_path_tmp, "popfit_final.RData", sep=""))
				popfit_final_old <- popfit_final
				load(file=paste(output_path_tmp, "popfit_quant.RData", sep=""))
				popfit_quant_old <- popfit_quant

				##	NOTE: There seems to be a bug combining forests for models
				##		based on different numbers of observations.  Therefore, to
				##		get this to work correctly, we're going to set the proximity
				##		and predicted attributes to NULL and 0 respectively before
				##		combining:
				popfit_final_combined$proximity <- NULL
        popfit_final_combined$predicted <- 0
        popfit_final_old$proximity <- NULL
        popfit_final_old$predicted <- 0
				
        popfit_quant_combined$predicted <- 0
        popfit_quant_old$predicted <- 0
				
				popfit_final_combined <- combine(popfit_final_combined, popfit_final_old)
        popfit_quant_combined <- combine(popfit_quant_combined, popfit_quant_old)
			}
			
			index <- index + 1
		}
		
		popfit_final <- popfit_final_combined
		popfit_quant <- popfit_quant_combined
		
		output_path_tmp <- paste(root_path, "/output/", country, "/tmp/", sep="")
		
		##  Save off our popfit object for this set of data:
		save(popfit_final, file=paste(output_path_tmp, "popfit_final_combined.RData", sep=""))
		#load(file=paste(output_path_tmp, "popfit_final_combined.RData", sep=""))
		save(popfit_quant, file=paste(output_path_tmp, "popfit_quant_combined.RData", sep=""))
		#load(file=paste(output_path_tmp, "popfit_quant_combined.RData", sep=""))
	}
}

##	END:	Estimate RandomForest model
#####



#####
##	BEGIN:	Predict for gridded covariates

##	Let's parallelize the process using the snow package:
##		NOTE:  The 00.1, alpha, script used the new parallel package
##			but it seems like this is a simpler way given that we have
##			to write out our results one block at a time...

##		NOTE also:  If you've changed the predictor set then you 
##			need to change the column renaming in the cluster_predict()
##			function and the subset of the RasterLayer objects in the
##			raster brick that gets created before the process is started...

##	Create the cluster process within a function:
cluster_predict <- function(prediction_raster, quant_output=FALSE, ...) {
	##	Start the timer:
	start_time = proc.time()[3]

	##	Pull the cluster:
	cl <- getCluster()
	on.exit( returnCluster() )
	
	##	Determine the number of cores we're working with:
	nodes <- length(cl)

	##	Generate a set of blocks on which to process the raster
	#blocks <- blockSize(prediction_raster, chunksize=200000, minblocks=nodes*4)
	blocks <- blockSize(prediction_raster, chunksize=100000, minblocks=nodes*4)

	pb <- tkProgressBar(title = "Predicting Population Density:", min = 0, max = blocks$n, width = 300)

	##	Pass off required libraries and data to the cluster workers:
	clusterEvalQ(cl, {
		require(raster)
		require(randomForest)
	})

	#clusterExport(cl, c("popfit", "covariate_stack", "transY"))
  if (quant_output) {
	  clusterExport(cl, c("popfit_final", "popfit_quant", "covariate_stack", "transY"))
  } else {
    clusterExport(cl, c("popfit_final", "covariate_stack", "transY"))
  }

	##	Since "blocks" only exists inside this function's environment, we need
	##		to call environment() to get our current environment rather than the 
	##		global environment assumed by clusterExport():
	clusterExport(cl, "blocks", envir=environment())
	clusterExport(cl, "quant_output", envir=environment())
	
	##	Define the function that will be run on each cluster to do the predictions:
	clFun <- function (i) {

		row_data <- data.frame( getValues(covariate_stack, row=blocks$row[i], nrows=blocks$nrows[i]) )
			
		##	Convert field names to something more manageable and that matches our popfit variable list:

		###	Full covariate stack:
		names(row_data) <- c(names(popfit_final$forest$xlevels), "census_mask", "water_raster")

		##	Convert to factors where appropriate:
		for (var_name in names(popfit_final$forest$xlevels)) {
			if (grepl("cls", var_name, perl=T)) {
				row_data[[var_name]] <- factor(row_data[[var_name]], levels=0:10, labels=as.character(0:10))
			}
		}

		##	Detect if we have any NA or Inf values, and that the values are 
		##		covered by our census administrative units:
		na_present <- apply(is.na(row_data), 1, any)
		inf_present <- apply(row_data == -Inf | row_data == Inf, 1, any)
		census_mask <- (is.na(row_data$census_mask))
		water_mask <- (row_data$water_raster == 1)

		##	Use the first if you want to mask out water pixels, this can greatly
		##		speed up predictions over areas with a lot of water, however, you
		##		run the risk of having no predictions in the resulting dataset
		##		if you have a census block small enough that it might only have
		##		water cover (GeoCover/GlobCover is what determines the water mask):
		roi_subset <- (!na_present & !inf_present & !census_mask & !water_mask)
		#roi_subset <- (!na_present & !inf_present & !census_mask)
			
		##	Create a set of predictions based on our covariates:
		predictions <- numeric(length=length(row_data[,1]))
		predictions[] <- NA

		#predictions <- data.frame("rf_pred"=predictions, "rf_pred_alt"=predictions, "rf_sd"=predictions, "rf_cv"=predictions)
		#predictions <- data.frame("rf_pred"=predictions, "rf_sd"=predictions)
		predictions <- data.frame("rf_pred"=predictions, "rf_sd"=predictions, "rf_05"=predictions, "rf_50"=predictions, "rf_95"=predictions)

		##	If we have data where NAs or Inf values are not present then we predict for those cells (where we subset our data according to the roi_subset and remove the census zone and water mask columns (length(row_data) - 2):
		if (sum(roi_subset) > 0) {
			#predictions[roi_subset] <- predict(popfit, newdata=row_data[roi_subset,1:(length(row_data)-2)])
			
			prediction_set <- predict(popfit_final, newdata=row_data[roi_subset,1:(length(row_data)-2)], predict.all=TRUE)
			
			predictions$rf_pred[roi_subset] <- transY(apply(prediction_set$individual, MARGIN=1, mean), inverse=TRUE)
			
			#predictions$rf_pred_alt[roi_subset] <- apply(transY(prediction_set$individual, inverse=TRUE), MARGIN=1, mean)
			#predictions$rf_sd[roi_subset] <- apply(transY(prediction_set$individual, inverse=TRUE), MARGIN=1, sd)
			#predictions$rf_cv[roi_subset] <- (predictions$rf_sd[roi_subset] / predictions$rf_pred[roi_subset]) * 100
			
			predictions$rf_sd[roi_subset] <- apply(prediction_set$individual, MARGIN=1, sd)
			
      if (quant_output) {
			  prediction_set <- predict(popfit_quant, newdata=row_data[roi_subset,1:(length(row_data)-2)], quantiles=c(0.05, 0.50, 0.95))
			  predictions$rf_05[roi_subset] <- transY(prediction_set[,1], inverse=TRUE)
			  predictions$rf_50[roi_subset] <- transY(prediction_set[,2], inverse=TRUE)
			  predictions$rf_95[roi_subset] <- transY(prediction_set[,3], inverse=TRUE)
		  }
		}

		return(predictions)
	}

	##	Start all nodes on a prediction:
	for (i in 1:nodes) {
		sendCall(cl[[i]], clFun, i, tag=i)
	}

	##	Start the raster writer object so we can store our results as they
	##		come back from our cluster:
	setwd(output_path)
	
	prediction_raster <- writeStart(prediction_raster, filename="predict_density.tif", format="GTiff", datatype="FLT4S", overwrite=TRUE, options=c("COMPRESS=LZW"))
	
	#prediction_raster_alt <- prediction_raster
	#prediction_raster_alt <- writeStart(prediction_raster_alt, filename="predict_density_alt.tif", format="GTiff", datatype="FLT4S", overwrite=TRUE, options=c("COMPRESS=LZW"))
	sd_raster <- prediction_raster
	sd_raster <- writeStart(sd_raster, filename="predict_density_sd.tif", format="GTiff", datatype="FLT4S", overwrite=TRUE, options=c("COMPRESS=LZW"))
	#cv_raster <- prediction_raster
	#cv_raster <- writeStart(cv_raster, filename="predict_density_cv.tif", format="GTiff", datatype="FLT4S", overwrite=TRUE, options=c("COMPRESS=LZW"))

  if (quant_output) {
  	prediction_raster_05 <- prediction_raster
  	prediction_raster_05 <- writeStart(prediction_raster_05, filename="predict_density_05.tif", format="GTiff", datatype="FLT4S", overwrite=TRUE, options=c("COMPRESS=LZW"))
  	prediction_raster_50 <- prediction_raster
  	prediction_raster_50 <- writeStart(prediction_raster_50, filename="predict_density_50.tif", format="GTiff", datatype="FLT4S", overwrite=TRUE, options=c("COMPRESS=LZW"))
  	prediction_raster_95 <- prediction_raster
  	prediction_raster_95 <- writeStart(prediction_raster_95, filename="predict_density_95.tif", format="GTiff", datatype="FLT4S", overwrite=TRUE, options=c("COMPRESS=LZW"))
  }

	##	Create our primary cluster processing loop, recalling that we already
	##		have clusters running:
	cat("Total blocks to process: ", blocks$n, "\n")
	for (i in 1:blocks$n) {
		##	Receive results from a node:
		predictions <- recvOneData(cl)

		##	Check if there was an error:
		if (!predictions$value$success) {
			stop("ERROR: Cluster barfed...\n\n", predictions)
		}

		##	Which block are we processing:
		block <- predictions$value$tag
		cat("Received block: ", block, "\n")
		flush.console()

		##	Now store our predictions in our prediction
		#prediction_raster <- writeValues(prediction_raster, transY(predictions$value$value, inverse=TRUE), blocks$row[block])
		
		prediction_raster <- writeValues(prediction_raster, predictions$value$value$rf_pred, blocks$row[block])
		#prediction_raster_alt <- writeValues(prediction_raster_alt, predictions$value$value$rf_pred_alt, blocks$row[block])
		sd_raster <- writeValues(sd_raster, predictions$value$value$rf_sd, blocks$row[block])
		#cv_raster <- writeValues(cv_raster, predictions$value$value$rf_cv, blocks$row[block])
		
    if (quant_output) {
  		prediction_raster_05 <- writeValues(prediction_raster_05, predictions$value$value$rf_05, blocks$row[block])
  		prediction_raster_50 <- writeValues(prediction_raster_50, predictions$value$value$rf_50, blocks$row[block])
  		prediction_raster_95 <- writeValues(prediction_raster_95, predictions$value$value$rf_95, blocks$row[block])
    }


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
	
	#prediction_raster_alt <- writeStop(prediction_raster_alt)
	sd_raster <- writeStop(sd_raster)
	#cv_raster <- writeStop(cv_raster)
	
  if (quant_output) {
	  prediction_raster_05 <- writeStop(prediction_raster_05)
	  prediction_raster_50 <- writeStop(prediction_raster_50)
	  prediction_raster_95 <- writeStop(prediction_raster_95)
  }

	close(pb)
	print(paste("Elapsed Prediction Time:", proc.time()[3] - start_time, "seconds"))
	flush.console()

	return(prediction_raster)
}


##	Create a raster stack of our cropped covariates and the zonal_raster file which will allow us to restrict processing to just the areas within the boundaries of our census data area (NOTE: This should be changed here to match the covariates used in the estimation of the model, as well as the renaming applied in the cluster predict function.  This will speed processing up slightly especially if used on a subset of the predictors):

##	Set the extent of processing:
my_extent <- extent(zonal_raster)
#my_extent <- extent(zonal_raster, 1000, 3000, 2000, 4000)
#my_extent <- extent(((extent(zonal_raster)@xmax - extent(zonal_raster)@xmin)/2 - 50000), ((extent(zonal_raster)@xmax - extent(zonal_raster)@xmin)/2 + 50000), ((extent(zonal_raster)@ymax - extent(zonal_raster)@ymin)/2 - 50000), ((extent(zonal_raster)@ymax - extent(zonal_raster)@ymin)/2 + 50000))
#plot(zonal_raster, ext=my_extent)


##	Create a zone-based raster that will be used as in/out of census block
##		mask:
#census_mask <- crop( zonal_raster, my_extent )

##	By using the zonal layer which was rasterized from the census data,
##		we have seams between regions processed through the RF method
##		separately.  Instead, let's create a mask based on our buffered data:
##	Read in our buffer shapefile:
census_buffer_path <- paste(project_path, "Census/Derived", sep="")
census_buffer <- readOGR(dsn=census_buffer_path, layer="census_buffer")
census_mask <- rasterize(census_buffer, zonal_raster)
census_mask <- crop( census_mask, my_extent )

##	Re-set the extent based on the sometimes slightly different extent
##		post-cropping:
my_extent <- extent(census_mask)


##	Confirm that our extent matches the minimum extent for all covariates.
##		Small inconsistencies can arrive as a product of the ArcGIS 
##		processing so this is a necessary step otherwise the stacking process
##		will fail in some cases:
for (var_name in names(popfit$forest$xlevels)) {
	assign("tmp_raster", raster( covariates[[var_name]]$path ))
	if (xmin(tmp_raster) > xmin(my_extent)) { my_extent@xmin <- xmin(tmp_raster) }
	if (xmax(tmp_raster) < xmax(my_extent)) { my_extent@xmax <- xmax(tmp_raster) }
	if (ymin(tmp_raster) > ymin(my_extent)) { my_extent@ymin <- ymin(tmp_raster) }
	if (ymax(tmp_raster) < ymax(my_extent)) { my_extent@ymax <- ymax(tmp_raster) }
}

##	Re-crop census_mask:
census_mask <- crop(census_mask, my_extent)


##	Create raster objects from each covariate raster:
##	NOTE: We are cropping each layer here because as I've discovered
##		with very large rasters if you try to crop the resulting covariate
##		stack the crop() function will fail...  This adds a small amount of
##		time onto the raster/crop process but seems to be more solid.
##	NOTE: There are some inconsistencies with the way that ArcGIS manages
##		to match the extent and snapping of rasters for output from projection
##		and so originally I was using a follow-up Con() call.  This, however,
##		was also inconsistent across different machines, so I removed it, which
##		means that it's very possible that we have rasters differing in numbers
##		of rows and columns by one pixel, especially among the large rasters
##		being cropped down (BioClim, Lights, HydroSHEDS, etc.).  But, the good
##		news is that the crop() call below will fix this as long as your 
##		raster files have identical pixel dimensions.  In the latest version 
##		of the raster package there's also an error that gets thrown:
##			Error : argument "rattypes" is missing, with no default
##		I can't figure out where this is coming from (no traceback() is 
##		generated), and the resulting cropped images seems to overlay perfectly
##		so we'll just go with it.
for (var_name in names(popfit$forest$xlevels)) {
	print(paste("Stacking: ", var_name, sep=""))
	flush.console()
	assign(var_name, crop(raster( covariates[[var_name]]$path ), my_extent))
	print(paste("ncell:", eval(parse(text=paste("ncell(", var_name, ")", sep="")))))
}

##	We need to ensure the water mask is included in covariate stack:
water_raster <- crop(raster( covariates[["lan_cls210"]]$path ), my_extent )


##	Create a raster object based on our census data zones file to store our predictions:
prediction_raster <- crop(raster(zonal_raster), my_extent)


##	Occasionally there is an error with cropping, a bug from the raster
##		package that results in matching dimensions but differing extents
##		and I don't know why as both ele and ele_slope are identical in
##		spatial structure.  If there's a problem with stacking then uncomment
##		the checking code and the setting of extents by hand for any layers
##		needed:

###	Check extents:
#for (var_name in c(names(popfit$forest$xlevels), "census_mask", "water_raster", "prediction_raster")) {
#	print(paste("Stacking: ", var_name, sep=""))
#	print(eval(parse(text=var_name)))
#	print(eval(parse(text=paste("extent(", var_name, ") == extent(census_mask)", sep=""))))
#	flush.console()
#}
#
#for (var_name in c(names(popfit$forest$xlevels), "census_mask", "water_raster", "prediction_raster")) {
#	print(paste("Fixing Extent: ", var_name, sep=""))
#	eval(parse(text=paste("extent(", var_name, ") <- extent(census_mask)", sep="")))
#	flush.console()
#}


##	Now stack all of our covariates and masks together:
covariate_stack <- eval(parse( text=
	paste(
		"stack(
			c(", 
			paste( c(names(popfit$forest$xlevels), "census_mask", "water_raster"), collapse=","),
			")
		)"
	)
))


###	Alternatively, could try a brick() instead of stack() for speed of 
###		extraction, however building the brick may take a very long time
###		as all raster covariates would need to be combined in memory...
#raster_list <- eval(parse( text=
#	paste(
#		"list(",
#			paste( c(names(popfit$forest$xlevels), "census_mask", "water_raster"), collapse=","),
#		")",
#	sep="")
#))
#
#covariate_stack <- brick(raster_list)


beginCluster(n=2)
prediction_raster <- cluster_predict(prediction_raster, quant_output=FALSE)
endCluster()

#Received block:  16614 
#[1] "Elapsed Prediction Time: 38116.63 seconds"

##	10.58 hours for a raster of 16614 x 79999 or 132895386 cells on a 
##		dual core laptop.  It	takes just under twice as long to run on a 
##		single core, with perhaps	a half hour of overhead tacked on.  More
##		complex randomForest models will take longer to run (up to 6.5 hours
##		dual core on the LUECI desktops).


##	Save the workspace:
#save.image(file=paste(output_path, "tmp", "00.2_image.RData", sep="/"))
#load(file=paste(output_path, "tmp", "00.2_image.RData", sep="/"))


##	Return our working directory to the source folder:
setwd(paste(root_path, "/src", sep=""))


##	END:	Predict for gridded covariates
#####
