## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~               Easy way to create an NC file            ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Creator: Javier
## type: Generic example
## Open library
library(ncdf4)


## Get x and y vectors (dimensions)
Longvector = seq(-180, 180, length = 50)
Latvector = seq(-90, 90, length = 30)
## Define data
dataset = list(1:18250)

## Define the dimensions
dimX = ncdim_def("Long", "degrees", Longvector)
dimY = ncdim_def("Lat", "degrees", Latvector)
dimT = ncdim_def("Time", "days", 1:18250)

## Define missing value
mv = -9999

## Define the data
var1d = ncvar_def( "var1d", "units", dimX, mv, prec="double")
var2d = ncvar_def( "var2d", "units", list(dimX,dimY), mv, prec="double")
var3d = ncvar_def( "var3d", "units", list(dimX,dimY,dimT), mv, prec="double")

## Create the NetCDF file
## If you want a NetCDF4 file, explicitly add force_v4=T
nc = nc_create("writevals.nc", list(var1d, var2d, var3d))

## Write data to the NetCDF file
ncvar_put(nc, var3d, dataset[[1]], start=c(1, 1, 1),
          count=c(1, 1, length(dataset[[1]])))

## Close your new file to finish writing
nc_close(nc)
