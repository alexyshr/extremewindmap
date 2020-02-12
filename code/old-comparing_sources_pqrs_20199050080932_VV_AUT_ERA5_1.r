#For study area (Colombia), load ERA5 information (netCDF file: outfile_nc4c_zip9.nc) 
#with wind variable "10fg: 10-m Wind-Gust since previous postprocesing".

library(ncdf4)
(ncname <- "outfile_nc4c_zip9")
(filename <- paste("data/", ncname, ".nc", sep = ""))
ncin <- nc_open(filename)
print (ncin)
lon <- ncvar_get(ncin, "longitude")
lon
nlon = dim(lon)
nlon
lat <- ncvar_get(ncin, "latitude")
lat
nlat = dim(lat)
nlat
ntime <-  dim(ncvar_get(ncin, "time"))
variablename <- "fg10"
fg10.units <- ncatt_get(ncin, variablename, "units")
fg10.units
lonlat.unstack <- expand.grid(lon=as.numeric(lon), lat=as.numeric(lat))

#With unstacked array of ERA5 coordinates, create a rectangular array of sf points covering 
#study area (isdcolpoints), and representing ERA5 cell centers, or ERA5 stations.
#Create point sf with lat and lon and value as cell index
isdcolpoints = st_as_sf(lonlat.unstack, coords=1:2, crs=st_crs(4326))
isdcolpoints$value = 1:(nlon*nlat)
str(isdcolpoints)
#plot(isdcolpoints)


#Convert ERA5 sf points to stars object, but:
#- each point at center of the cell
#- as cell value an unique index (cell index), corresponding to the row index of lonlat.unstack
#Cell shape for this netCDF file is a square: 0.25º x 0.25º

library(stars)
pointsbbox = st_bbox(isdcolpoints)
cellsize = lonlat.unstack$lon[2]- lonlat.unstack$lon[1]

mybbox = st_bbox(c(pointsbbox$xmin - (cellsize/2), pointsbbox$xmax + (cellsize/2), pointsbbox$ymax + (cellsize/2), pointsbbox$ymin - (cellsize/2)), crs = st_crs(4326))
(isdcolraster.st = st_rasterize(isdcolpoints, st_as_stars(mybbox, nx = nlon, ny = nlat, values = isdcolpoints$value)))


#Plot ERA5 raster, ERA5 points, and stations sample IDEAM
library(RColorBrewer)
n <- 60 #Number of colors to take from palette
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#Plot raster and points
plot(isdcolraster.st, nbreaks=nlon*nlat, col=sample(col_vector, (nlon*nlat)-1, replace = TRUE), reset= FALSE, xlim=c(min(lon)-1,max(lon)+1), ylim=c(min(lat)-1,max(lat)+1))
plot(isdcolpoints, pch=".", cex=.5, col='black', add=TRUE)

#Plot with sample stations (IDEAM and ISD) and ERA5 raster.
plot(isdcolraster.st, nbreaks=nlon*nlat, col=sample(col_vector, (nlon*nlat)-1, replace = TRUE), reset= FALSE, xlim=c(-80,-66), ylim=c(-6,14))
plot(st_geometry(stationssample_isd), pch=16, cex=1, col='black', add=TRUE)
plot(st_geometry(stationssample_ideam), pch=21, cex=0.5, col='white', bg='black',add=TRUE)
text(st_coordinates(stationssample_ideam), labels = paste(1:length(stationssample[,1]), stationssample_ideam$codigo1, sep='-'), pos = 2, cex=0.8)


