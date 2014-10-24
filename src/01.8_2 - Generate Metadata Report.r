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
project_path <- paste(root_path, "data/", country, "/", sep="")
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
##	BEGIN:	Package loading

require(knitr)
require(markdown)

##	BEGIN:	Package loading
#####



#####
##	BEGIN:	Report generation
start_time = proc.time()[3]

##	Change the working directory to the /src folder to load in a list of 
##		all R markdown (.Rmd) files:
src_path <- paste(root_path, "src/", sep="")
report_files <- list.files(src_path, pattern='\\.Rmd$', full.names = FALSE)

##	Change the working directory back to the output temporary folder for
##		the creation of our knitted output:
setwd(output_path_tmp)
for (report_file in report_files) { 
  knit2html( paste(src_path, report_file, sep="") ) 

	knitted_report_name <- paste(output_path_tmp, sub("\\.Rmd", ".html", report_file), sep="")
	output_report_name <- paste(output_path, country, "_metadata.html", sep="")
	if (file.exists(output_report_name)) {
		unlink(output_report_name)
	}
	file.copy(knitted_report_name, output_report_name) 
}
print(paste("Elapsed Processing Time:", proc.time()[3] - start_time, "seconds"))


##	If you encounter an error in processing or need to fix the .md file
##		by hand (located in the output/tmp folder), then you can re-run
##		the knitting process of just the .md file after fixing:

#if (file.exists(output_report_name)) {
#	unlink(output_report_name)
#}
#markdown::markdownToHTML(file=paste(output_path_tmp, sub("\\.Rmd", ".md", report_file), sep=""), output=output_report_name)


##	Return our working directory to the source folder:
setwd(paste(root_path, "src", sep=""))


##	END:	Report generation
#####
