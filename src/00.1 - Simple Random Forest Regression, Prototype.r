require(parallel)
require(foreign)
require(randomForest)
require(yaImpute)
require(raster)
##	NOTE: The tcltk library is used to create a progress bar, but the
##		raster package has a progress bar mechanism built in in recent versions.
##		Look up the pbCreate(), pbStep(), and pbClose() functions in the raster
##		package...
require(tcltk)


setwd("D:/Documents/Graduate School/Research/Population/Data/Modeling Work/data")

census_data <- read.dbf("Afghanistan/AFG-popdata_w_lc.dbf")


summary(census_data)


census_data$pop_dens <- census_data$EST0910 / census_data$HECTARES

y_data <- census_data$pop_dens
x_data <- census_data[,21:29]

##	Because our summary from the "Tabulate" tool in ArcGIS is the sum of 
##		the area in decimal degrees for each land cover class inside each
##		polygon, we need to perform an additional calculation to generate
##		percentage or proportion areas.  Luckily there is a Shape_Area field
##		that should be the shape area in decimal degrees as well...  We need
##		to divide each land class sum by that total area:
for (i in 1:9) {
	x_data[,i] <- x_data[,i] / census_data$Shape_Area
}

##	Confirm that we have proportions that sum to 1.0 for each row (each
##		polygon within the census data):
sum(x_data[1,])
summary(apply(x_data, 1, sum))



##########
##	Random Forest regression of population density:
set.seed(1002)

##	Check our memory footprint just to make sure:
memory.size()


##	Now we are going tune our randomForest population density regression: 
tuneRF(x=x_data, y=y_data, plot=TRUE, mtryStart=4, ntreeTry=1000, improve=0.0001, stepFactor=1.5, trace=TRUE) 

##	mtry = 4  OOB error = 35.30821 
##	Searching left ...
##	mtry = 3        OOB error = 41.33215 
##	-0.1706100 1e-04 
##	Searching right ...
##	mtry = 6        OOB error = 29.58690 
##	0.1620393 1e-04 
##	mtry = 9        OOB error = 22.43372 
##	0.2417683 1e-04 
##	  mtry OOBError
##	3    3 41.33215
##	4    4 35.30821
##	6    6 29.58690
##	9    9 22.43372


#	Based on the above I'll use an mtry value of 3:
start_time = proc.time()[3]
popfit = randomForest(x=x_data, y=y_data, mtry=9, nodesize=5, ntree = 1000, na.action = na.omit, importance = TRUE, keep.forest=TRUE)
print(paste("Elapsed Fitting Time:", proc.time()[3] - start_time, "seconds"))
popfit

##	Call:
##	 randomForest(x = x_data, y = y_data, ntree = 1000, mtry = 9,      nodesize = 5, importance = TRUE, keep.forest = TRUE, na.action = na.omit) 
##	               Type of random forest: regression
##	                     Number of trees: 1000
##	No. of variables tried at each split: 9
##	
##	          Mean of squared residuals: 20.96704
##	                    % Var explained: 67.6

summary(popfit)
plot(popfit)


##	Check our diagnostics:
##	Recall that running predict on the randomForest object will return OOB predictions:
predict(popfit)
importance(popfit)

##	Get the variable names of the top 20 predictors:
names(importance(popfit)[,"%IncMSE"][order(importance(popfit)[,"%IncMSE"], decreasing=TRUE)])[1:20]

varImpPlot(popfit)
varUsed(popfit)
summary(treesize(popfit))
partialPlot(x=popfit, pred.data=x_data, x.var="VALUE_190")
plot(popfit)


##	For continuous regression, plot observed vs. predicted:
plot(x=y_data, y=predict(popfit))
abline(a=0, b=1, lty=2)

##	For continuous regression, plot residuals vs. observed:
plot(x=y_data, y=(y_data - predict(popfit)))
abline(a=0, b=0, lty=2)

##	For continuous regression, plot residuals vs. fitted:
plot(x=predict(popfit), y=(y_data - predict(popfit)))
abline(a=0, b=0, lty=2)
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



##########
##	Create density estimates on 30m grid from land cover:

setwd("D:/Documents/Graduate School/Research/Population/Data/MDA Land Cover/Refined")

land_cover <- raster("afghanistan_lc_mosaic_reclass_missing_filled_rurb_8bit.img")

#plot(land_cover)
cells <- cellFromRowCol(land_cover, 1, 1)
cells <- cellFromRowColCombine(land_cover, 17000:17005, 18000:18005)
sample_data <- land_cover[cells]
xyFromCell(land_cover, cells)


predict_density <- function(values, na.rm=FALSE) {
	##	In our case, 0 values are NAs, so we convert them:
	values[values == 0] <- NA

	##	If our na.rm flag is TRUE remove any NAs from the value list:
	if (na.rm) values <- values[!is.na(values)]
	
	##	Check to see if we have NAs in the value list and if so return NA:
	if (sum(is.na(values)) > 0) {
		return(NA)
	}

	##	Generate counts of each land cover type:
	counts <- table(values)

	##	Generate a total sum of cells:
	total <- sum(counts)

	##	Create a data.frame object to hold our "predictor" values:
	new_x_data <- data.frame( 
		"VALUE_11"=0,
		"VALUE_40"=0,
		"VALUE_130"=0,
		"VALUE_140"=0,
		"VALUE_160"=0,
		"VALUE_190"=0,
		"VALUE_200"=0,
		"VALUE_210"=0,
		"VALUE_240"=0
	)

	##	Generate the proportion of each land cover within the focal window
	##		and store it in the data.frame:
	for (i in c("11", "40", "130", "140", "160", "190", "200", "210", "240")) {
		if(!is.na(counts[i])) new_x_data[[paste("VALUE_", i, sep="")]][1] <- counts[i]/total
	}

	##	Predict using our random_forest or yai prediction using the new_x_data
	##		and return it as the focal cell value:

	##	yai impute (note rownames must be unique and new for new data for
	##		yai imputation to succeed):
	#rownames(new_x_data) = paste(rownames(new_x_data),"_new",sep="")
	#prediction <- impute(newtargets(popfit, newdata=new_x_data))$y

	##	randomForest prediction:
	prediction <- predict(popfit, newdata=new_x_data)


	return(prediction)
}

##	Run our predict function on our 5x5 sample:
predict_density(sample_data)


##	Now use our focal function to predict population density at a 30m
##		scale across our image.  We are using a 5x5 neighborhood matrix
##		to calculate the focal window:
setwd("D:/Documents/Graduate School/Research/Population/Data/Modeling Work/output")

##	Set a subset for testing:
land_cover_subset <- crop(land_cover, extent(c(70.25, 70.50, 34.50, 34.75)))
#land_cover_subset <- crop(land_cover, extent(c(70.45, 70.50, 34.70, 34.75)))
plot(land_cover_subset)



##	NOTE: So this approach will take waaaaay too long, so we're going to 
##		create a grid at roughly 100 m, sample our cells at each grid location
##		and then generate our predictions at each point:
#start_time = proc.time()[3]
##prediction_image <- focal(land_cover_subset, w=5, fun=predict_density, pad=FALSE, na.rm=FALSE, filename="prediction_image.tif", overwrite=TRUE, datatype="FLT4S", progress="window")
#prediction_image <- focal(land_cover_subset, w=5, fun=predict_density, pad=FALSE, na.rm=FALSE, progress="window")
#print(paste("Elapsed Prediction Time:", proc.time()[3] - start_time, "seconds"))
#
#summary(prediction_image)
#plot(prediction_image)



##	Create a coordinate set that we are going to loop over to extract cell
##	values for and create estimates of population density centered on:

##	Can't load the whole image in as a SpatialPixelsDataFrame... too big!
#land_cover_sp <- readGDAL("afghanistan_lc_mosaic_reclass_missing_filled_rurb_8bit.img")

##	NOTE: This process is hugely speeded up by tiling it, extracting a chunk
##		of the data into memory at a time and processing the whole thing from
##		memory rather than accessing cell blocks from the on-disk raster:

land_cover_extent <- extent(land_cover_subset)

x_list <- seq(land_cover_extent@xmin + res(land_cover_subset)[1]*100/30/2, land_cover_extent@xmax - res(land_cover_subset)[1]*100/30/2, by=res(land_cover_subset)[1]*100/30)
y_list <- seq(land_cover_extent@ymin + res(land_cover_subset)[2]*100/30/2, land_cover_extent@ymax - res(land_cover_subset)[2]*100/30/2, by=res(land_cover_subset)[2]*100/30)

predictions <- raster(nrows=length(y_list), ncols=length(x_list), ext=land_cover_extent, crs=projection(land_cover_subset))



###	Create our coordinate list and an empty vector to hold our predictions:
coords <- expand.grid(x=x_list, y=y_list)

##	For testing purposes of the parallelized function below:
i <- 1

##	TODO: So the creation of a prediction raster that has the same extent as
##		the land cover subset (tile) that you're working with doesn't really work
##		as the corner coordinates using a 5x5 window of values will spill over into
##		NA territory (focal_x and focal_y for the corner coordinate is only 1 pixel
##		away from the edge, not two...  So for the larger loop creation the subset
##		tile that is read from disk into memory for the land_cover file should
##		have a border created...


##	Create a function that's built to predict for a set of coordinates:
predict_value <- function(i) {
	focal_x <- colFromX(land_cover_subset, coords[i,1])
	focal_y <- rowFromY(land_cover_subset, coords[i,2])
	cell_list <- cellFromRowColCombine(land_cover_subset, rownr=(focal_y-2):(focal_y+2), colnr=(focal_x-2):(focal_x+2) )
	
	sample_data <- land_cover_subset[cell_list] 

	return(predict_density(sample_data))
}

##	Create a function that can operate over a list of coordinate indices
##		which is what the clusterization algorithm will want:
predict_at_coordinates <- function(indices) {
	##	Create progress bar:
	pb <- tkProgressBar(title = "Predicting Population Density:", min = min(indices), max = max(indices), width = 300)

	##	For the progress bar:
	min_index <- min(indices)
	max_index <- max(indices)

	##	Run the loop over our coordinates and create our predictions:
	##		NOTE: So, because we are going to clusterize this process, we need
	##			to have our output function kick out only a single object...
	prediction_values <- data.frame("prediction"=NA)
	prediction_values[length(coords[,1]),1] <- NA
	for (i in indices) {
		prediction_values[i,1] <- predict_value(i)
		setTkProgressBar(pb, i, label=paste( round((i - min_index)/(max_index - min_index)*100, 0),"% done"))
	}
	close(pb)
	return(prediction_values)
}



##	Test the above function on a single core:
start_time = proc.time()[3]

prediction_values <- predict_at_coordinates(1:length(coords[,1]))

##	And now create a raster object from our predictions and coordinates:
prediction_raster <- SpatialPixelsDataFrame(points=SpatialPoints(cbind(coords[,1], coords[,2])), data=data.frame("prediction"=prediction_values))
prediction_raster <- raster(prediction_raster)
start_time = proc.time()[3]

prediction_values <- predict_at_coordinates(1:length(coords[,1]))

##	And now create a raster object from our predictions and coordinates:
prediction_raster <- SpatialPixelsDataFrame(points=SpatialPoints(cbind(coords[,1], coords[,2])), data=data.frame("prediction"=prediction_values))
prediction_raster <- raster(prediction_raster)
#plot(prediction_raster)

print(paste("Elapsed Prediction Time:", proc.time()[3] - start_time, "seconds"))



##	Let's parallelize the process using the new parallel package:
cl <- makeCluster(getOption("cl.cores", 2))
#cl <- makeCluster(getOption("cl.cores", 1))

##	Create the sequence list to run by core:
cluster_sequences <- clusterSplit(cl=cl, seq=1:length(coords[,1]))

##	Set up the cluster environment:
clusterEvalQ(cl, {
	require(raster)
	require(randomForest)
	require(tcltk)
})

clusterExport(cl, c("popfit", "predict_density", "land_cover_subset", "coords", "predict_value", "predict_at_coordinates"))

##	Test cluster output:
clusterEvalQ(cl, {
	predict_value(10000)
})

##	Send the coordinate list to our two cores:
start_time = proc.time()[3]
cluster_output <- clusterApply(cl=cl, fun="predict_at_coordinates", x=cluster_sequences)

##	Merge the pieces:
prediction_values <- cluster_output[[1]]
for(i in 2:length(cluster_output)) {
	prediction_values[is.na(prediction_values)] <- cluster_output[[i]][is.na(prediction_values)]
}
prediction_raster <- SpatialPixelsDataFrame(points=SpatialPoints(cbind(coords[,1], coords[,2])), data=data.frame("prediction"=prediction_values))
prediction_raster <- raster(prediction_raster)
plot(prediction_raster)


print(paste("Elapsed Prediction Time:", proc.time()[3] - start_time, "seconds"))

##	Stop the cluster:
stopCluster(cl)



##	This technique should be roughly the equivalent to the above, using the 
##		coordinates to extract out cells based on a particular point location,
##		and indeed seems just about as fast...

#land_cover_tile <- as.matrix(land_cover_subset)
#predictions <- land_cover_tile
#predictions[,] <- NA
#
###	Create progress bar:
#start_time = proc.time()[3]
#pb <- tkProgressBar(title = "Predicting Population Density:", min = 0, max = length(predictions), width = 300)
#
#count <- 0
#
###	Run the loop over our coordinates and create our predictions:
###	Process every 3rd focal cell:
#for (i in seq(3, (nrow(land_cover_tile)-2), by=3)) {
#	for (j in seq(3, (ncol(land_cover_tile)-2), by=3)) {
##for (i in 3:(nrow(land_cover_tile)-2)) {
##	for (j in 3:(ncol(land_cover_tile)-2)) {
#		sample_data <- land_cover_tile[(i-2):(i+2), (j-2):(j+2)]
#		predictions[i,j] <- predict_density(sample_data)
#		
#		count <- count + 3	
#		setTkProgressBar(pb, count, label=paste( round(count/length(predictions)*100, 0),"% done"))
#	}
#}
#close(pb)
#print(paste("Elapsed Prediction Time:", proc.time()[3] - start_time, "seconds"))
