REM Set these to reflect the location where you have unzipped the RF working
REM 	folders and the location and version of Python/ArcGIS for your particular
REM 	environment...

set RFHOME=D:\Documents\Graduate School\Research\Population\Analysis\RF\RF.1.0
set ARCPY=C:\Python26\ArcGIS10.0


REM In practice, nothing below this should need to be set...


set PATH=%PATH%;%ARCPY%\;%ARCPY%\Scripts;%RFHOME%\bin\R\R-3.0.1\bin\x64;%RFHOME%\bin\FWTools2.4.7\bin;%RFHOME%\bin\SAGA-GIS;%RFHOME%\bin\SAGA-GIS\modules;

set R_LIBS_USER=%RFHOME%\bin\R\win-library\3.0
set R_HOME=%RFHOME%\bin\R\R-3.0.1\
set R_PATH=%RFHOME%\bin\R\R-3.0.1\


.\RStudio\bin\rstudio.exe ..\src\