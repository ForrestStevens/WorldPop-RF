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
##	NOTICE:  In practice nothing below this line should need to be regularly
##		edited unless very specific details need to be changed about the 
##		modeling process (e.g. speciyfing an extent or something detailed
##		about the RandomForest modeling).
#####



#####
##	BEGIN:	Package loading and configuration

require(raster)
require(GSIF)
require(plotKML)
require(rgdal)
require(R.utils)

##	Occasionally RSAGA won't load, so you may need to require it on 
##		its own:
require(RSAGA)


##	Set up data paths:


##	First, we check to see if we used an alternative census file for map
##		production than the one we used for RF model parameterization:
##	TODO:  Eventually this needs to be pulled from the Metadata.R file
##		instead of a hard coded directory presence.
dataset_folders <- list.dirs(path=project_path, full.names=FALSE, recursive=FALSE)

if (paste(project_path, "/! New Census", sep="") %in% dataset_folders) {
	census_folder_name <- "! New Census"
	new_census_section <- TRUE
} else {
	census_folder_name <- "Census"
	new_census_section <- FALSE
}



dataset_name <- list.files(paste(project_path, census_folder_name, "/", sep=""),  "shp$")
census_data <- readOGR(dsn=paste(project_path, census_folder_name, sep=""), substr(dataset_name, 1, nchar(dataset_name)-4))
census_year <- census_data[["YEARPOP"]][1]


##	Default coarse level output (for automated output from my scripts):
output_file <- paste(output_path, country, "_ppp_v", rf_version, "_", census_year, ".tif", sep="")


kml_output_path <- paste(output_path_tmp, "kml/", sep="")
dir.create(kml_output_path, showWarnings=FALSE)


##	END:	Package loading and configuration
#####



#####
##	BEGIN:	Internal function declaration


##	We are going to overload the png() function to be the resolution of our 
##		underlying data:
png <- function(width=480, height=480, ...) {
	if (width==480 & height==480) {
		##	Since we're now using our tiled processing, we are going to kick
		##		out 100x100 pixel tiles so hard code our width and height to match
		##		the 100 * 0.0008333 degree tiles:
		#if ((raster_data@ncols * raster_data@nrows) > (3000 * 3000)) {
		#	scale_factor <- sqrt((3000^2) / (raster_data@ncols * raster_data@nrows))
		#	width <- round(raster_data@ncols*scale_factor)
		#	height <- round(raster_data@nrows*scale_factor)
		#} else {
		#	width <- raster_data@ncols
		#	height <- raster_data@nrows
		#}
		width = 500
		height = 500
	}

	add_args <- modifyList(list(type="cairo-png"), list(...))

	##	Always use the cairo-png driver so that we always use 24-bit palette:
	return(grDevices::png(width=width, height=height, ...))
}


##	This is my own version of plotKML with a signature of "list", edited to
##		handle a list of RasterLayer objects:
frs_plotKML <- function (obj, ...) {
	.local <- function (obj, folder.name = normalizeFilename(deparse(substitute(obj, env = parent.frame()))), file.name = paste(folder.name, ".kml", sep = ""), size = NULL, colour, points_names = "", shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png", plot.labpt = TRUE, labels = "", metadata = NULL, kmz = get("kmz", envir = plotKML.opts), ...) {
		if (any(!(sapply(obj, class) == "SpatialPointsDataFrame" | sapply(obj, class) == "SpatialLinesDataFrame" | sapply(obj, class) == "SpatialPolygonsDataFrame" | sapply(obj, class) == "SpatialPixelsDataFrame" | sapply(obj, class) == "RasterLayer"))) {
			stop("List of objects of class SpatialPoints*, SpatialLines*, SpatialPolygons*, SpatialPixels*, RasterLayer expected")
		}
		kml_open(folder.name = folder.name, file.name = file.name)
		
		for (i in 1:length(obj)) {
			if (missing(colour)) {
				##	Catch for RasterLayer objects to use the name of the raster layer:
				if (class(obj[[i]]) != "RasterLayer") {
					obj[[i]]@data[, "colour"] <- obj[[i]]@data[, 1]
				}
			} else {
				if (is.name(colour) | is.call(colour)) {
					obj[[i]]@data[, "colour"] <- eval(colour, obj[[i]]@data)
				} else {
					obj[[i]]@data[, "colour"] <- obj[[i]]@data[, as.character(colour)]
				}
			}
		}
		
		if (all(sapply(obj, class) == "SpatialPointsDataFrame")) {
			for (i in 1:length(obj)) {
				if (points_names == "") {
					if (is.numeric(obj[[i]]@data[, 1])) {
						points_names_i <- signif(obj[[i]]@data[, 1], 3)
					}	else {
						points_names_i <- paste(obj[[i]]@data[, 1])
					}
				}

				if (is.numeric(obj[[i]]@data[, "colour"])) {
					kml_layer.SpatialPoints(obj[[i]], colour = colour, points_names = points_names_i, shape = shape, metadata = metadata, ...)
				} else {
					kml_layer.SpatialPoints(obj[[i]], colour = colour, points_names = points_names, shape = shape, metadata = metadata, ...)
				}
			}
		}

		if (all(sapply(obj, class) == "SpatialLinesDataFrame")) {
			for (i in 1:length(obj)) {
				kml_layer.SpatialLines(obj[[i]], metadata = metadata, ...)
			}
		}

		if (all(sapply(obj, class) == "SpatialPolygonsDataFrame")) {
			for (i in 1:length(obj)) {
				if (labels == "") {
					obj[[i]]@data[, "labels_i"] <- obj[[i]]@data[, 1]
				} else {
					if (is.name(labels) | is.call(labels)) {
						obj[[i]]@data[, "labels_i"] <- eval(labels, obj[[i]]@data)
					} else {
						obj[[i]]@data[, "labels_i"] <- obj[[i]]@data[, deparse(labels)]
					}
				}
			
				kml_layer.SpatialPolygons(obj[[i]], colour = colour, plot.labpt = plot.labpt, labels = labels_i, metadata = metadata, ...)
			}
		}

		if (all(sapply(obj, class) == "SpatialPixelsDataFrame")) {
			for (i in 1:length(obj)) {
				bbn <- round(diff(obj[[i]]@bbox[1, ])/obj[[i]]@grid@cellsize[1]) * round(diff(obj[[i]]@bbox[2, ])/obj[[i]]@grid@cellsize[2])
				if (max(obj[[i]]@grid.index, na.rm = TRUE) > bbn) {
					x <- as.data.frame(obj[[i]])
					suppressWarnings(gridded(x) <- ~x + y)
					proj4string(x) = obj[[i]]@proj4string
					obj[[i]] <- x
				}

				raster_name_i <- paste(names(obj[[i]])[1], "_", i, ".png", sep = "")
				kml_layer.SpatialPixels(obj[[i]], colour = colour, raster_name = raster_name_i, metadata = metadata, plot.legend = FALSE, ...)
			}
		}


		##	Create a new class within the list processing script to handle RasterLayers:
		if (all(sapply(obj, class) == "RasterLayer")) {
			plot_legend = TRUE
			for (i in 1:length(obj)) {
				##	Suppressing any warnings from data with all NAs:
				if (suppressWarnings(cellStats(obj[[i]], stat="max")) > 0) {
					##	Substitute any zeroes with NA:
					#obj[[i]] <- subs(obj[[i]], data.frame("from"=0, "to"=NA), subsWithNA=FALSE)
					##	Convert current RasterLayer object:
					raster_name_i <- paste(names(obj[[i]]), "_", i, ".png", sep = "")
					kml_layer.Raster(obj[[i]], raster_name = raster_name_i, metadata = metadata, plot.legend = plot_legend, ...)
					plot_legend <- FALSE
				}
			}
		}


		kml_close(file.name = file.name)
		if (kmz == TRUE) {
			kml_compress(file.name = file.name)
		}
		
		kml_View(file.name)
	}
	
	.local(obj, ...)
}


##	We are also going to override the internal function from GSIF that tries
##		to use the Windows registry to locate FWTools so that we can use our
##		own, hardcoded path to the FWTools folder.  This uses the handy
##		reassignInPackage() function from the R.utils package:
frs_FWTools.path <- function(...) { return(paste(root_path, "/bin/FWTools2.4.7/", sep="")) }
reassignInPackage(".FWTools.path", pkgName="GSIF", frs_FWTools.path)


##	Last, we are going to actually define our own tile function instead of
##		using the one that exists in the GSIF package, because we want to 
##		control the way that the tile image is saved (the GSIF version uses
##		gdalwarp and loses the NA values, converting them to zeroes:
#frs_tile <- function (x, y, block.x, tmp.file = FALSE, program, show.output.on.console = FALSE, ...) {
#	if (filename(x) == "") {
#		stop("Function applicable only to 'RasterLayer' objects linking to external raster files")
#	}
#
#	if (missing(program)) {
#		if (.Platform$OS.type == "windows") {
#			fw.dir = .FWTools.path()
#			program = shQuote(shortPathName(normalizePath(file.path(fw.dir, "bin/gdalwarp.exe"))))
#		} else {
#			program = "gdalwarp"
#		}
#	}
#
#	if (missing(y)) {
#		b <- bbox(x)
#		pol <- SpatialPolygons(list(Polygons(list(Polygon(matrix(c(b[1, 1], b[1, 1], b[1, 2], b[1, 2], b[1, 1], b[2, 1], b[2, 2], b[2, 2], b[2, 1], b[2, 1]), ncol = 2))), ID = "1")), proj4string = CRS(proj4string(x)))
#		y <- getSpatialTiles(pol, block.x = block.x, return.SpatialPolygons = FALSE, ...)
#	}
#
#	x.lst <- list()
#	message("Clipping raster object using 'gdalwarp'...")
#
#	for (j in 1:nrow(y)) {
#		if (tmp.file == TRUE) {
#			outname <- tempfile()
#		} else {
#			outname <- paste(normalizeFilename(deparse(substitute(x, env = parent.frame()))), j, sep = "_")
#		}
#
#		##	TODO: We need to set the -srcnodata and -dstnodata flags to use the
#		##		-Inf (-2E38) value for the Raster object being fed into it.  This 
#		##		tile function is not any different than the standard GSIF function
#		##		for	Raster data at the moment:
#		try(system(paste(program, shortPathName(normalizePath(filename(x))), set.file.extension(outname, ".tif"), "-te", y[j, 1], y[j, 2], y[j, 3], y[j, 4]), show.output.on.console = show.output.on.console))
#		try(x.lst[[j]] <- raster(set.file.extension(outname, ".tif")))
#	}
#
#	return(x.lst)
#}

##	NOTE: Instead of using the above tile function, modified from the GSIF
##		package, I'm going to use the below function instead.  This function
##		relies on the Raster package crop() function, and keeps cropped 
##		raster objects in memory rather than writing them out to disk.  This 
##		has the advantage of being faster and not triggering the problem of
##		NODATA values being written out to the output by gdalwarp as zeroes
##		by default.  If we run into a country where this in-memory processing
##		is too problematic, then we can modify the TODO section above to use
##		gdalwarp, but enforce proper handling of NODATA values:
frs_raster_tile <- function (x, y, block.x, tmp.file = FALSE, program, show.output.on.console = FALSE, ...) {
	if (filename(x) == "") {
		stop("Function applicable only to 'RasterLayer' objects linking to external raster files")
	}

	if (missing(y)) {
		b <- bbox(x)
		pol <- SpatialPolygons(list(Polygons(list(Polygon(matrix(c(b[1, 1], b[1, 1], b[1, 2], b[1, 2], b[1, 1], b[2, 1], b[2, 2], b[2, 2], b[2, 1], b[2, 1]), ncol = 2))), ID = "1")), proj4string = CRS(proj4string(x)))
		y <- getSpatialTiles(pol, block.x = block.x, return.SpatialPolygons = FALSE, ...)
	}

	x.lst <- list()
	message("Clipping raster object using raster...")

	for (j in 1:nrow(y)) {
		x.lst[[j]] <- crop(x, extent(y[j,1],y[j,3],y[j,2],y[j,4]))
	}

	return(x.lst)
}


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
#setwd(paste(project_path, census_folder_name, "/Derived", sep=""))

#census_data <- readOGR(".", "census_covariates")
raster_data <- raster(output_file)

##	END:	Data loading and data structure setup
#####



#####
##	BEGIN:	Plot an figure for KML conversion


##	Sample plot using a custom palette:

#plot(raster_data, zlim=c(0,35), col=colorRampPalette(c("navyblue", "steelblue", "limegreen", "yellow", "white"))(255))


setwd(kml_output_path)
data(SAGA_pal)


##	Set the zlim values for legend and plotting:
z.lim = c(0,35)

##	This will fail on large countries with big rasters... Therefore we need
##		to tile them using the following:
#plotKML(raster_data, file.name=paste(country, "_ppp_v", rf_version, "_", census_year, ".kml", sep=""), colour_scale=SAGA_pal[[1]], z.lim=z.lim)

##	We no longer need to do this because we are tiling our output at 100 
##		pixel wide/high, but this is what you need to do if you regenerate
##		the PNG being mapped to Google Earth:

###	Now re-plot our image at full resolution:
#png(filename = paste(country, "_ppp_v", rf_version, "_", census_year, ".kml", sep=""), bg = "transparent", type = "cairo-png")
#par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
#obj <- raster_data
#obj <- calc(obj, fun = function(x) {
#	x[x < z.lim[1]] <- z.lim[1]
#	return(x)
#})
#obj <- calc(obj, fun = function(x) {
#	x[x > z.lim[2]] <- z.lim[2]
#	return(x)
#})
#
##raster::image(obj, col=SAGA_pal[[1]], zlim=z.lim, frame.plot=FALSE, main="", maxpixels=raster_data@ncols*raster_data@nrows)
##raster::image(obj, col=SAGA_pal[[1]], zlim=z.lim, frame.plot=FALSE, main="", maxpixels=3000^2)
#raster::image(obj, col=colorRampPalette(c("navyblue", "steelblue", "limegreen", "yellow", "white"))(255), zlim=z.lim, frame.plot=FALSE, main="", maxpixels=3000^2)
#dev.off()


##	We are going to create a set of tiles for our country, that have a 
##		maximum of 500 * 0.0008333 degrees width/height so that we can be
##		sure to have maximum resolution output for our Google Earth KML.  Note
##		that the block.x parameter is in the units of the projection, in our
##		case degrees:
setwd(kml_output_path)
#raster_data_list <- tile(raster_data, block.x=500 * 0.0008333)
raster_data_list <- frs_raster_tile(raster_data, block.x=500 * 0.0008333)


#showMethods("plotKML")
#getMethod("plotKML", signature="list")

##	Use our custom function with RasterLayer added to the list processing
##		options:
##	Define color scale for the full list, then replace the first one with
##		#FFFFFF which is the transparent color so that zeroes are clear:
color_vec <- colorRampPalette(c("navyblue", "steelblue", "limegreen", "yellow", "#FEFEFE"))(255)
#color_vec[1] <- "#FFFFFF"
frs_plotKML(raster_data_list, file.name=paste(country, "_ppp_v", rf_version, "_", census_year, ".kml", sep=""), colour_scale=color_vec, z.lim=z.lim, kmz=FALSE)


##	Now the last thing that we need to do is to convert all of our output 
##		PNG files to 8-bit depth so that Google Earth doesn't drop a load 
##		due to transparency issues (4-bit and 2-bit PNG with transparent 
##		colors don't work in Google Earth).  To do this we use a system
##		call to ImageMagick's mogrify command:
shell(paste('"', root_path, '\\bin\\ImageMagick-6.7.4-Q16\\mogrify.exe" -format png -depth 8 -type TruecolorMatte -define png:color-type=6 *.png', sep=""), intern=TRUE)
shell(paste('"', root_path, "\\bin\\7-Zip\\7z.exe\" -Tzip a ..\\..\\", country, "_ppp_v", rf_version, "_", census_year, ".kmz *.*", sep=""), intern=TRUE)



##	Example:
#data(eberg_grid)
#gridded(eberg_grid) <- ~x+y
#proj4string(eberg_grid) <- CRS("+init=epsg:31467")
#data(SAGA_pal)
#plotKML(eberg_grid["TWISRT6"], colour_scale = SAGA_pal[[1]])


##	Return our working directory to the source folder:
setwd(paste(root_path, "/src", sep=""))


##	END:	Plot an figure for KML conversion
#####

