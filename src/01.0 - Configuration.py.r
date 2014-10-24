##	This is the primary per-configuration file, parsed by each script,
##		in this folder, including both R and Python scripts.  There should
##		be fewer per-script configuration items now as these will control
##		the country-name and root path (path to the folder containing the RF
##		/src, /data, /output, etc. folders.

##	NOTE: Make sure that you don't use the R assignment operator <- as
##		this file is also parsed by Python, so use = instead:


##	Configure the country abbreviation and name:
country = "IDN"

##	The version of the scripts used to produce the mapping products, and
##		which will match the "_v" portion of the filename outputs:
rf_version = "2b"

##	This should be set to the folder *containing* the "RF" folder structure,
##		but make sure that you use forward slashes instead of back slashes or
##		double back slashes:
root_path = "D:/Documents/Graduate School/Research/Population/Analysis/RF/Working RF/"
