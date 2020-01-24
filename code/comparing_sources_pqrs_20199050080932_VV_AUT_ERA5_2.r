#Get datetime object from netCdf file 
t.units <- ncatt_get(ncin, "time", "units")
t.units
time.array <-  ncvar_get(ncin, "time")
nt = dim(time.array)
nt
head(time.array, 1)
tail(time.array, 1)
tustr <- strsplit(t.units$value, " ")
#tustr
#Our data is in Gregorian Calendar: "hours"      "since"      "1900-01-01" "00:00:00.0"
anio_mes_dia = unlist(tustr)[3]
anio_mes_dia
#Note that the dates are converted first to seconds!
library(lubridate)
timestamp = as_datetime(c(time.array*60*60),origin=anio_mes_dia)
timestamp_string <- as.character(timestamp)
str(timestamp_string)

#Get slices from netCDF

for(i in 1:length(stationssample[,1])){
  #for(i in 1:1){  
  assign(paste("isd_fg10", i, isd_ideam_intersect["ideam_id",i], sep = "_"), 
         ncvar_get(ncin, 'fg10', start=c(as.integer(isd_ideam_intersect["ideam_lon_index_era5",i]),as.integer(isd_ideam_intersect["ideam_lat_index_era5",i]),1), count=c(1,1,ntime)))    
}

#Create xts objects (ERA5)
library(xts)
#station_slice79_12.5 = as.xts(station_slice79_12.5)
for(i in 1:length(stationssample[,1])){
  #for(i in 1:1){
  assign(paste("station_slice", i, isd_ideam_intersect["ideam_id",i], sep = "_"), xts(x= eval(parse(text = (paste("isd_fg10", i, isd_ideam_intersect["ideam_id",i], sep = "_")))), order.by = timestamp))    
}