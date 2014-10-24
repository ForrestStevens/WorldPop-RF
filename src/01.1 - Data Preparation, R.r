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
country <- "NGA"

##	This should be set to the folder *containing* the "RF" folder structure:
root_path <- "D:/Documents/Graduate School/Research/Population/Data/"
project_path <- paste(root_path, "RF/data/", country, "/", sep="")

##	Load the metadata from the file created in the country's /data folder:
source(paste(project_path, "Metadata.r", sep=""))


##	It's good practice to check the folders in your country's data folder
##		to make sure they match the names specified in your Metadata.r file.
##		Compare what's printed to those folders in the list.  If any folders
##		are not present in data folder, but specified in the metadata names
##		make sure that they are "DEFAULT" datasets that we construct on the
##		fly (see the "default_datasets" configuration parameter below).  Also
##		be sure that any folders in the data folder have a corresponding
##		entry in the Metadata.r list, otherwise they will not be processed
##		in the Random Forest estimation.  (NOTE:  You can exclude any folder 
##		or file from further processing by placing an exclamation point in the
##		file name...):
names(metadata)[order(names(metadata))]

##	Specify required datasets that *must* have data provided in the
##		folder of the same name.  This should *never* change but it is 
##		included in the configuration options here as a reminder to 
##		double check their presence in both the country's data folder *and* an
##		entry in the Metadata.r file.  As of now, there must be a shapefile, 
##		named anything you like in the Census folder, with ADMINID and ADMINPOP
##		fields, and	an *.img or *.tif file in the Landcover that has at least 
##		a subset of the land covers as coded in our paper:
required_datasets <- c("Census", "Landcover")

##	Specify the datasets with default data (specialized treatment using
##		global datasets and hand-processed from mosaics or existing data
##		that we will use.  In practice these won't change much, if ever,
##		however, occasionally you may want to exclude one or more, like NPP,
##		from default processing.  Because they are included here *does not*
##		mean that you can't override the default processing and dataset
##		options by dropping your own ShapeFile or raster dataset in the folder.
##		But you may exclude default processing altogether by removing the name
##		from this list *or* by removing the name from the Metadata.r file in 
##		the country's /data folder:
default_datasets <- c("NPP", "Lights", "Roads", "Rivers", "Waterbodies", "Populated", "Protected", "Elevation", "Temp", "Precip" )


##	END:	Load configuration options
#####



#####
##	NOTICE:  In practice nothing below this line should need to be regularly
##		edited unless very specific details need to be changed about the 
##		modeling process (e.g. speciyfing an extent or something detailed
##		about the RandomForest modeling).
#####



#####
##	BEGIN: Package loading and fixed parameter configuration

##	This is the R JSON package for exporting our metadata and covariate
##		configuration data into XML for Python to import:
require(rjson)


##	The other packages require for data pre-processing are for automated
##		MODIS data download:

#detach(package:MODIS, unload=TRUE)
#remove.packages("MODIS")

##	From R-Forge:
#install.packages("MODIS", repos="http://R-Forge.R-project.org", type="source")
##	Try this one first:
#install.packages("MODIS", repos="http://R-Forge.R-project.org")

##	Alternative, install from local:
#install.packages("D:/Programming/R/MODIS Processing in R/MODIS Package SVN/pkg/MODIS", repos=NULL, type="source")
#install.packages("D:/Programming/R/MODIS Processing in R/MODIS Package SVN, Backup After My Edits/pkg/MODIS", repos=NULL, type="source")

##	For package updated and created for source:
#setwd(source_path)
#install.packages("../lib/MODIS.tar.gz", repos=NULL, type="source")

require(MODIS)
require(rgdal)


##	Check to make sure we have all dependencies:
MODIS:::checkDeps()

##	> setRepositories()
##	--- Please select repositories for use in this session ---
##	
##	
##	1: + CRAN
##	2: + CRAN (extras)
##	3:   BioC software
##	4:   BioC annotation
##	5:   BioC experiment
##	6:   BioC extra
##	7: + Omegahat
##	8: + R-Forge
##	9:   rforge.net
##	
##	Enter one or more numbers separated by spaces, or an empty line to cancel
##	1: 1 2 7 8
##	
##	install.packages(c("bitops", "rgeos", "maps", "mapdata", "maptools", "ptw", "XML", "SSOAP", "XMLSchema"), dependencies=TRUE)


##	Double check the dependencies are good by running he checkDeps() again
##		and potentially installing some from source...


##	Last, double check that our tool chain is set up correctly for the MODIS
##		Reprojection Toolkit and the GDAL version:
MODIS:::checkTools("GDAL")
#MODIS:::checkTools("MRT")


##	MODIS local archive output location, as specified in your .MODIS_Opts.R 
##		file using the outDirPath configuration line.  This file should
##		be located in your documents folder in Windows, e.g.:
##			C:\Users\Forrest\Documents\.MODIS_Opts.R file:
modis_out <- paste(options()$MODIS_outDirPath, "/", country, sep="")


##	END: Package loading and fixed parameter configuration
#####



#####
##	BEGIN:	Internal function declaration

##	END:	Internal function declaration
#####



#####
##	BEGIN:	Metadata setup and compilation

##	Next create a list containing information about each of our derived raster datasets used as covariates in our models:
covariates <- list()
var_names <- character()


##	Get a list of folders that don't begin with an exclamation point, and
##		assume that they are covariate folders.  We'll keep track of them
##		and any that are left over after processing the "Metadata.r" file
##		we'll warn the user about, because the Python script will try to 
##		process any	of those folders:
dataset_folders <- dir(path = project_path, pattern = "^[^!]", full.names = FALSE, recursive = FALSE)
dataset_folders <- dataset_folders[grepl("^[^(Metadata.r)]", dataset_folders, perl=TRUE)]


##	Recall that metadata is loaded in from the "Metadata.r" file contained
##		in our per-country /data folder:
for (dataset in metadata) {

	##	Perform a variety of checks of the metadata and data before continuing:


	##	Remove the dataset folder from the list of folders in our directory:
	dataset_folders <- dataset_folders[dataset$dataset_folder != dataset_folders]


	##	Default variable name, folder location and dataset name for the 
	##		datasets:
	var_name <- substr( tolower(dataset$dataset_folder), 1, 3)
	dataset_base_path <- paste(project_path, dataset$dataset_folder, "/", sep="")
	dataset_derived_path <- paste(dataset_base_path, "Derived/", sep="")
	dataset_name <- tolower(dataset$dataset_folder)

	##	Add our three-letter variable name extracted from the folder name
	##		to a list and check to make sure there are no duplicates.  If there
	##		are throw an error and immediately bail:
	if (var_name %in% var_names) {
		stop(paste("ERROR:  It appears as though there is already a variable set and hence folder name that begins with \"", var_name, "\" and these must be unique.  Please rename the ", dataset$dataset_folder, " folder and update the Metadata.r file appropriately and re-run this script before proceeding!", sep=""))
	} else {
		var_names <- append(var_names, var_name)
	}

	##	As we iterate through our metadata datasets, check to make sure
	##		the folder exists and if it is not part of the default data
	##		to process there is a file (ShapeFile, *.img, *.tif) corresponding
	##		to the data type in the folder:
	dir.create(dataset_base_path, showWarnings=FALSE)
	dir.create(dataset_derived_path, showWarnings=FALSE)

	
	##	If the dataset is non-default or required then it *must* have a data
	##		file in the root dataset_folder:
	if (dataset$dataset_class == "raster") {
		dataset_name = list.files(dataset_base_path, "tif$|img$")[1]

		if (is.na(dataset_name) & ((dataset$dataset_folder %in% required_datasets) | !(dataset$dataset_folder %in% default_datasets))) {
			stop(paste("ERROR:  The ", dataset$dataset_folder, "\" is either a required raster dataset or has no global default and requires that you place a dataset in the folder prior to further processing.  Please update the Metadata.r file appropriately and re-run this script before proceeding!", sep=""))
		}
	} else {
		dataset_name = list.files(dataset_base_path, "shp$")[1]

		if (is.na(dataset_name) & ((dataset$dataset_folder %in% required_datasets) | !(dataset$dataset_folder %in% default_datasets))) {
			stop(paste("ERROR:  The ", dataset$dataset_folder, "\" is either a required linear or polygon dataset or has no global default and requires that you place a dataset in the folder prior to further processing.  Please update the Metadata.r file appropriately and re-run this script before proceeding!", sep=""))
		}
	}

	
	##	The last check we make is to ensure that the file dataset_name 
	##		specified in the Metadata.r file matches the file name we pulled
	##		above.  Just because it exists in a non-default dataset does not
	##		mean it is recorded and specified correctly.  This is important 
	##		because we want the metadata reports to be accurate:
	if (!is.na(dataset_name)) {
		if (dataset_name != dataset$dataset_name) {
			stop(paste("ERROR:  The ", dataset$dataset_name, " specified for the \"", dataset$dataset_folder, "\" does not match the file found in the folder on the drive:  ", dataset_name, ".  Please double check that the Metadata.r file has been updated correctly to match the data processing you really want in order to ensure that the metadata reports are correct.", sep=""))
		}
	}
		
		
	##	If we've passed all checks above for the current dataset then proceed
	##		to build the derived covariates:

	derived_sets <- substr( dataset$derived, 1, 3)

	for (derived in dataset$derived[ derived_sets %in% c("cls", "dst", "prp", "", "slo") ]) {
		var_name <- substr( tolower(dataset$dataset_folder), 1, 3)
		dataset_name <- tolower(dataset$dataset_folder)
		
	
		##	No longer needed as this was simplified in the newer metadata and
		##		Python processing script:

		###	Build up dataset_name variable from options and derived products:
		#if ("merged" %in% dataset$derived) {
		#	dataset_name <- paste( tolower(dataset_name), "_merged", sep="")
		#}


		if (derived != "") {
			dataset_name <- paste( dataset_name, "_", derived, ".tif", sep="")
		} else {
			dataset_name <- paste( dataset_name, ".tif", sep="")
		}
		
		switch( substr(derived, 1,3),
			cls = {
				var_name <- paste( var_name, "_", derived, sep="")
				dataset_description = paste("Binary classification of", tolower(dataset$dataset_folder), "as a factor")
				dataset_summary = "modal"
			},
			dst = {
				var_name <- paste( var_name, "_", derived, sep="")
				dataset_description = paste("Distance to feature or value of", tolower(dataset$dataset_folder))
				dataset_summary = "mean"
			},
			prp = {
				var_name <- paste( var_name, "_", derived, sep="")
				dataset_description = paste("Proportion of grid cells within polygon or five cell radius of grid location that are", tolower(dataset$dataset_folder))
				dataset_summary = "mean"
			},
			slo = {
				var_name <- paste( var_name, "_", derived, sep="")
				dataset_description = "Estimated slope from DEM"
				dataset_summary = "mean"
			},
			{
				##	Note: that for non-summarized datasets we just set the output to match the name of the folder without the derived output (it may be blank):
				dataset_description = paste("Mean value of", tolower(dataset$dataset_folder), "within polygon or value at grid location")
				dataset_summary = "mean"
			}
		)

		covariates[[var_name]] <- list(
			dataset_folder      = dataset$dataset_folder,
			dataset_name        = dataset_name,
			dataset_description = dataset_description,
			dataset_summary     = dataset_summary
)
		covariates[[var_name]][["path"]] <- paste(dataset_derived_path, covariates[[var_name]]$dataset_name, sep="")
	}
}


##	Save these to structures on disk as they are used in the metadata report
##		scripts:
save(metadata, file=paste(output_path_tmp, "metadata.RData", sep=""))
save(covariates, file=paste(output_path_tmp, "covariates.RData", sep=""))

##	We are also going to export these metadata lists to JSON file objects
##		so that Python can load them:
metadata_json <- toJSON(metadata)
covariates_json <- toJSON(covariates)

cat(metadata_json, file=paste(output_path_tmp, "metadata.json", sep=""))
cat(covariates_json, file=paste(output_path_tmp, "covariates.json", sep=""))

if (length(dataset_folders) > 0) {
	message <- paste("WARNING:  The following folders and/or unknown files were detected in the country data folder without a proper entry in the processed \"Metadata.r\" file.  Please double check them as they will be processed by the \"Data Preparation, Python.py\" file but won't be included in the Random Forest model without a proper Meatadata.r entry.  If you do not intend for them to be processed then please add an exclamation point to the beginning of their name or move them outside the country's data folder:\n", paste(dataset_folders, collapse=", "))
	warning(message)
}

##	END:	Metadata setup and compilation
#####



#####
##	BEGIN:	Download and process MODIS products:


##	If NPP is specified in the metadata as a dataset to include download
##		and process it to create our output:
if (("NPP" %in% names(metadata)) & ("NPP" %in% default_datasets)) {

	##	Load data associated with the extent of our current country:

	setwd(paste(project_path, "Census/", sep=""))
	census_file_name <- list.files(".", "shp$")

	census <- readOGR(".", layer=substr(census_file_name,1,nchar(census_file_name)-4))


	##	Double-check the requirement that our Census file is in Geographic
	##		coordinates (because we have to hardcode the pixel size for 
	##		reprojection to ~100m:
	if (proj4string(census) != "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") {
		stop("ERROR:  The census data shapefile must be projected to GCS WGS84 in order to proceed!")
	}


	##	To just download the tiles:
	#getHdf(product="MOD17A3", begin="2010-01-01", end="2011-01-01", extent=extent(census), localArcPath=modis_path)

	##	To check the layer names of a downloaded file:
	#getSds(HdfName="D:/tmp/MODIS/MOD17A3.055/2010.01.01/MOD17A3.A2010001.h21v08.055.2011276120715.hdf")

	
	##	Set the out_SDSstring to by just for NPP, the middle of the three layers:
	out_SDSstring <- "010"


	##	To download (if not downloaded already, the local archive is checked) and process tiles using GDAL we can do the following:

	##	NOTE: runGdal() requires that we have FWTools or some other software
	##		with a working, binary implementation of GDAL that R can find:
	runGdal(job=country, product="MOD17A3", begin="2010-01-01", end="2011-01-01", extent=extent(census), SDSstring=out_SDSstring, outProj="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", pixelSize=0.0008333, resamplingType="NN")

	##	NOTE: There's a bug in the runGdal() function that tries move a file
	##		after it is complete, which throws a "Warning message:" like the 
	##		following.  It is safe to ignore but make sure to check the NPP
	##		folder to make sure the MODIS output look correct.
	##	Warning message:
	##	In file.remove(ofile) :
 	##	 cannot remove file 'C:\tmp\MODIS_ARC\PROCESSED\KHM\MOD17A3.A2010001.Npp_1km.tif', reason 'No such file or directory'


	##	Since we had to send output to a temporary directory you can now move the contents of the temporary output folder into the /data folder for the country in question.
	dataset_name <- list.files(modis_out, "tif$|img$", full.names=FALSE)
	file.copy(from = paste(modis_out, dataset_name, sep="/"), to = paste(project_path, "NPP/", dataset_name, sep=""), overwrite=TRUE)

}


##	END:	Download and process MODIS products:
#####
