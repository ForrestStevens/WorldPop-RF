## @knitr DataLoad

require(raster)
require(randomForest)
##	NOTE: The tcltk library is used to create a progress bar, but the
##		raster package has a progress bar mechanism built in in recent versions.
##		Look up the pbCreate(), pbStep(), and pbClose() functions in the raster
##		package...
require(tcltk)


load("D:/Documents/Graduate School/Research/Population/Data/RF - Prototype/output/AFG.old/00.2_image.RData")


##	Only required if running by hand and not read in from a knitr report:
if (!exists("path")) path <- "D:/Documents/Graduate School/Research/Population/Data/RF/"
if (!exists("country_prefix")) country_prefix <- "KHM"

print(path)
print(country_prefix)

##	This is unneccessary because it should be sourced from the output 
##		folder during report creation and a relative path should work:
#setwd(paste(path, "data", country_prefix, sep="/"))

##	This is a comment...
#x = rnorm(100)
#boxplot(x)
#plot(x,x)
#paste("The P-value is", t.test(x)$p.value)



print(popfit)
importance(popfit)[order(importance(popfit)[,1], decreasing=TRUE),]
varImpPlot(popfit)

###	For continuous regression, plot predicted vs. observed:
plot(x=y_data, y=predict(popfit), ylim=c(0,120), xlim=c(0,120))
abline(a=0, b=1, lty=2)

###	For continuous regression, plot residuals vs. observed:
plot(x=y_data, y=(y_data - predict(popfit)), xlim=c(0,20))
abline(a=0, b=0, lty=2)
