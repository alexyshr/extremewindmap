Sys.setenv(TZ='UTC')

#Stations equivalence:
#stationssample = data.frame(
#  isd_usaf_id= as.character(
#    c(803980, 803700, 802110, 802100, 801120, 801100, 800970, 800940, 800630, 800360, 800350, 800280)),
#  ideam_id = as.character(
#    c(48015050, 52055230, 26125061, 26125710, 23085270, 27015330, 16015501, 23195502, 13035501, 28025502, 15065180, 29045190)),
#  fc_zo =
#    c(1.197230052, 1.219771719, 1.102205474 , 1.18154867, 1.113341504, 1.29241596, 1.102205474, 1.177562503, 1.114968586, 1.184849704, 1.5, 1.224744381),
#  stringsAsFactors=FALSE)

stationssample = data.frame(
  isd_usaf_id= as.character(
    c(800360, 800350)),
  ideam_id = as.character(
    c(28025502, 15065180)),
  fc_zo =
    c(1.184849704, 1.5),
  stringsAsFactors=FALSE)


#fc_zo =
#  c(1.197230052, 1.219771719, 1.102205474 , 1.18154867, 1.113341504, 1.29241596, 1.102205474, 1.177562503, 1.114968586, 1.184849704, 1.143077647, 1.224744381),

#GMT
#Sys.timezone()

#vv_aut_2 in mts/seg
#
path_vv_aut_2 = "../data/"

#path_vvag_media_d = "E:/Thesis/ideamdata/PQRS_20199050100232/"

#myFile = read.table(paste0(path_vvag_media_d, "VVAG_MEDIA_D@11025010.data"), header=TRUE, sep="|", stringsAsFactors=FALSE) #this is for vvag_media_d
#myFile = read.table(paste0(path_vv_aut_2, "VV_AUT_2@13035501.data"), header=TRUE, sep="|", stringsAsFactors=FALSE)
#str(myFile)

#filenames = list.files(path_vv_aut_2, pattern = "VVAG_MEDIA_D@") #for vvag_media_d
filenames = list.files(path_vv_aut_2, pattern = "VV_AUT_2@") #for VV_AUT_2
#r = lapply(filenames, read.table, header=TRUE, sep="|", stringsAsFactors=FALSE)

ldf = lapply(filenames, function(x) {
  #dat = read.table(paste0(path_vvag_media_d, x), header=TRUE, sep="|", stringsAsFactors=FALSE) #this is for vvag_media_d
  dat = read.table(paste0(path_vv_aut_2, x), header=TRUE, sep="|", stringsAsFactors=FALSE)

  # Add column names
  names(dat) = c('fecha', 'valor')

  # Add a column with the year
  #dat$station_id = substring(x,14,nchar(x)-5)  #this is for vvag_media_d
  #Take control of 999.9, 314.00, 99.900, 1000
  dat$valor[dat$valor > 300] = NA
  dat$valor[dat$valor < 1] = NA

  dat$station_id = substring(x,10,nchar(x)-5)  #this is for vv_aut_2
  dat$mydatetime = as.POSIXct(dat$fecha,format="%Y-%m-%d %H:%M:%S", tz="UTC")
  return(dat)
})  #be aware that ldf is a list of dataframes and you need to do [[]], to go inside it

str(ldf)

#lapply(ldf, summary)

#Calculate Hourly Mean of each dataframe (station) inside ldf, then create ldf_hourlymean
#And unstack from rows to columns
ldf_hourlymean <- NULL
for(station in ldf){
  library(xts)
  library(dplyr)
  select <- dplyr::select
  myxts = na.omit(xts(x=select(station, valor), order.by = station$mydatetime))
  #Conversion to km/h (times 3.6) - Here we have V3600 (Hourly Mean)
  #myxts$valor = myxts$valor * 3.6

  #Conversion to 3-s gust using Durst curve
  #Gust factor from Durst curve = 1.52222
  myxts$valor = myxts$valor * 1.52


  #Correction factor by Roughness
  fczo = stationssample$fc_zo[stationssample$ideam_id == station$station_id[1]]
  #for (v in 1:length(stationssample$isd_usaf_id)) {
  #  if(station$station_id[1] == stationssample[v, 2]){
  #    fczo = stationssample[v, 3]
  #    break
  #  }
  #}
  myxts$valor = myxts$valor * fczo

  #temps_hours <- split(myxts, f = "hours")
  #temps_avg <- lapply(X = temps_hours, FUN = mean)

  colnames(myxts) = paste0("X", station$station_id[1])
  #myxts$station_id = station$station_id
  #myxts2 = to.period(myxts, period="hours") #Not meaning
  endhour = endpoints(myxts,on="hours")
  myxtshour = xts::period.apply(myxts, INDEX=endhour, FUN=mean)
  index(myxtshour)=trunc(index(myxtshour),"hours")   ###DOWN DOWN DOWN
  #myxtshour = align.time(myxtshour, n=3600) # UP UP UP #Write only YYYY-MM-DD HH:00:00, rounding to the next hour
  ldf_hourlymean <- cbind(ldf_hourlymean, myxtshour)
}


#myxtshour2 = align.time(myxtshour - lubridate::minutes(60), n=60*60)

#df = do.call("rbind", ldf) #with do.call is not any more a list of dataframe (only one dataframe)
#str(df)

#Conversion to km/h (times 3.6) - Here we have V3600 (Hourly Mean)
#ldf_hourlymean$valor = ldf_hourlymean$valor * 3.6

#Conversion to 3-s gust using Durst curve
#Gust factor from Durst curve = 1.52222
#ldf_hourlymean$valor = ldf_hourlymean$valor * 1.52


#
library('RPostgreSQL')
pg = dbDriver("PostgreSQL")
con1 = dbConnect(pg, user="user1", password="user1", host="localhost", port=5432, dbname="winddata")

#Load ISD Stations from database, but filtering to retrieve only the sample stations
wherestring3 = stationssample$isd_usaf_id
wherestring3 = paste("'", wherestring3, "'", sep = "")
wherestring3 = paste(wherestring3, collapse = ", " , sep = " ")
wherestring3 = paste("usaf IN (", wherestring3, ")", sep ="")

originalfields3 = c("id", "usaf", "latitud", "longitud")
originalfields3 = paste(originalfields3, collapse= ", ", sep = "")
query3 = paste("select", originalfields3, "from isd_all_stations", "where", wherestring3, sep=" ")
isd_stations = tbl(con1, sql(query3))

#Load IDEAM Stations from database, but filtering to retrieve only the sample stations
wherestring4 = stationssample$ideam_id
wherestring4 = paste("'", wherestring4, "'", sep = "")
wherestring4 = paste(wherestring4, collapse = ", " , sep = " ")
wherestring4 = paste("codigo1 IN (", wherestring4, ")", sep ="")

originalfields4 = c("objectid", "codigo1", "latitud", "longitud")
originalfields4 = paste(originalfields4, collapse= ", ", sep = "")
query4 = paste("select", originalfields4, "from ideam_all_stations", "where", wherestring4, sep=" ")
ideam_stations = tbl(con1, sql(query4))


#ISD stations and IDEAM stations to 'simple feature' (sf) object with points. Only spatial points that match station sample list.
library(sf)
isd_stations = st_as_sf(as_tibble(isd_stations), coords = c("longitud", "latitud"), crs = 4326)
ideam_stations = st_as_sf(as_tibble(ideam_stations), coords = c("longitud", "latitud"), crs = 4326)
stationssample_isd = isd_stations
stationssample_ideam = ideam_stations

#ERA5 stuff
#
#
#
#
#
#
source('comparing_sources_pqrs_20199050080932_VV_AUT_ERA5_1.r')



# With the sample stations points (ISD, IDEAM), get the raster index of ERA5 file. Raster index can be:
#   - row index (for longitud):
#   ISD: isd_lon_index
# IDEAM: ideam_lon_index
# - col index (for lat):
#   ISD: isd_lat_index
# IDEAM: ideam_lat_index
# - longitud
# ISD: isd_lon
# IDEAM: ideam_lon
# - latitud
# ISD: isd_lat
# IDEAM: ideam_lat
# - Cell index (corresponding to the row index of lonlat.unstack)
# ISD: isd_station
# IDEAM: ideam_station
#
# All this values are stored by rows in matrix isd_ideam_intersect. Columns represents the sample points.
#
# Other rows in this matrix:
# - isd_intersect_order: Order of the intersect operation between ISD sf points and raster. Depends on order of points in sf. Points in sf are not in the same order comparing ISD and IDEAM.
# - ideam_intersect_order: Order of the intersect operation between IDEAM sf points and raster. Depends on order of points in sf. Points in sf are not in the same order comparing ISD and IDEAM.
# - ideam_id: To test the match operation done to get raster indexes, ideam_id is calculated twice (using ISD points, and IDEAM points). Both values must be the same.
# - usaf: To test the match operation done to get raster indexes, usaf is calculated twice (using ISD points, and IDEAM points). Both values must be the same.

#Index due to ISD stations
isdcellindex = st_intersects(isdcolraster.st,stationssample_isd)
isdcellindexmatrix = as.matrix(isdcellindex)
isdlista1 = which(isdcellindexmatrix)
(isdrowlista1 = isdlista1%%length(lonlat.unstack[, 1]))
(isdlatindex = ceiling(isdrowlista1/nlon))
(isdlonindex = (isdrowlista1%%nlon))

isd_intersect = rbind(isd_lon_index_era5 = isdlonindex, isd_lat_index_era5 = isdlatindex, isd_station_era5 = isdrowlista1)
isd_intersect = rbind(usaf = as.character(stationssample_isd$usaf), isd_intersect )
isd_match = match(isd_intersect["usaf", ], stationssample$isd_usaf_id) #Find indexes in second parameter
isd_intersect = rbind(ideam_id = (as.character(stationssample$ideam_id[isd_match])),isd_intersect)
isd_intersect = rbind(isd_intersect_order = 1:length(stationssample[,1]), isd_intersect)
isd_intersect = isd_intersect[, order(as.integer(isd_intersect["ideam_id",]), decreasing = FALSE)]

ideamcellindex = st_intersects(isdcolraster.st,stationssample_ideam)
ideamcellindexmatrix = as.matrix(ideamcellindex)
ideamlista1 = which(ideamcellindexmatrix)
(ideamrowlista1 = ideamlista1%%length(lonlat.unstack[, 1]))
(ideamlatindex = ceiling(ideamrowlista1/nlon))
(ideamlonindex = (ideamrowlista1%%nlon))

ideam_intersect = rbind(ideam_lon_index_era5 = ideamlonindex, ideam_lat_index_era5 = ideamlatindex, ideam_station_era5 = ideamrowlista1)
ideam_intersect = rbind(ideam_id = as.character(stationssample_ideam$codigo1), ideam_intersect )
ideam_match = match(ideam_intersect["ideam_id", ], stationssample$ideam_id) #Find indexes in second parameter
ideam_intersect = rbind(usaf = as.integer((as.character(stationssample$isd_usaf_id[ideam_match]))),ideam_intersect)
ideam_intersect = rbind(ideam_intersect_order = 1:length(stationssample[,1]), ideam_intersect)
ideam_intersect = ideam_intersect[, order(as.integer(ideam_intersect["ideam_id",]), decreasing = FALSE)]

isd_lonlat=matrix(data=NA, nrow=2, ncol=length(stationssample[,1]))
for(i in 1:length(stationssample[,1])){
  isd_lonlat[1,i] = lonlat.unstack$lon[as.integer(isd_intersect["isd_station_era5",i])]
  isd_lonlat[2,i] = lonlat.unstack$lat[as.integer(isd_intersect["isd_station_era5",i])]
}
isd_intersect = rbind(isd_lon = as.character(isd_lonlat[1,]), isd_lat = as.character(isd_lonlat[2,]), isd_intersect)

ideam_lonlat=matrix(data=NA, nrow=2, ncol=length(stationssample[,1]))
for(i in 1:length(stationssample[,1])){
  ideam_lonlat[1,i] = lonlat.unstack$lon[as.integer(ideam_intersect["ideam_station_era5",i])]
  ideam_lonlat[2,i] = lonlat.unstack$lat[as.integer(ideam_intersect["ideam_station_era5",i])]
}

ideam_intersect = rbind(ideam_lon = as.character(ideam_lonlat[1,]), ideam_lat = as.character(ideam_lonlat[2,]), ideam_intersect)

(isd_ideam_intersect = rbind(isd_intersect, ideam_intersect))

#write.table(isd_ideam_intersect, file="isd_ideam_intersect13.csv", sep=",")


#Get ISD wind data (table isd_lite_unstack) from PostgreSql
#Important: Dataset is filtered using whereclause, retrieving only sample stations.
originalfields1 = stationssample$isd_usaf_id
#correction factor for Isd!!
fczo = stationssample$fc_zo
newfields1 = paste ("X", originalfields1, sep="")
originalfields1 = paste('"', originalfields1, '"', sep = "")
newfields1 = paste('"', newfields1, '"', sep = "")
fiedls_query1 = paste("CASE WHEN", originalfields1, "< 1", "THEN NULL ELSE", originalfields1, "*", fczo, "END AS", newfields1, sep = " ")

#No correction factors!
#fiedls_query1 = paste(originalfields1, "as", newfields1, sep = " ")
#Applying correction factors!!
#fiedls_query1 = paste(originalfields1, "*", fczo, "as", newfields1, sep = " ")

fiedls_query1 = c(paste('"', "mydatetime", '"', sep = ""), fiedls_query1)
fiedls_query1 = paste (fiedls_query1, "", sep= "", collapse=", ")
wherestring1 = stationssample$isd_usaf_id
wherestring1 = paste('"', wherestring1, '"', sep = "")

#Leaving small values
#wherestring1 = paste(wherestring1, "IS NOT NULL", sep = " ")
#wherestring1 = paste(wherestring1, collapse = " OR " , sep = " ")

#Not leaving small values (removing < 1) >>>Not valid, use CASE
wherestring1 = paste(wherestring1, ">= 1 AND", wherestring1, "IS NOT NULL", sep = " ")
wherestring1 = paste(wherestring1, collapse = ") OR (" , sep = " ")
wherestring1 = paste("(", wherestring1, ")", sep ="")

query1 = paste("select", fiedls_query1, "from isd_lite_unstack", "where", wherestring1, sep=" ")

cat(query1)
isdlite = tbl(con1, sql(query1))

#Datetime object for ISD
#???timestamp_isdlite <- as.POSIXct(as_tibble(select(isdlite, mydatetime))$mydatetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")
timestamp_isdlite <- as.POSIXct(as_tibble(select(isdlite, mydatetime))$mydatetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")


#str(timestamp_isdlite)
#head(timestamp_isdlite, 5)
#tail(timestamp_isdlite, 5)
#length(timestamp_isdlite)

#Ideam stations are ready


all_vvmx_aut_60 = as.data.frame(ldf_hourlymean) #Note that the variable name is related to all_vvmx_aut_60
                                                #but the data is related to vv_aut_2  >>>>>><<<<<< EYES!!!!!

all_vvmx_aut_60$mydatetime = rownames(all_vvmx_aut_60)
all_vvmx_aut_60 = select(all_vvmx_aut_60, mydatetime, starts_with("X"))
timestamp_all_vvmx_aut_60 <- as.POSIXct(as_tibble(select(all_vvmx_aut_60, mydatetime))$mydatetime,format="%Y-%m-%d %H:%M:%S", tz="UTC")

#Sort stations in isdlite and all_vvmx_aut_60 according to isd_ideam_intersect order.
#Although wind data from ISD and IDEAM (isdlite and all_vvmx_aut_60) are already filtered with sample stations, thus the stations in isdlite and all_vvmx_aut_60 are the same compared to the result of point-raster intersection array (isd_ideam_intersect).
#The selection only will change the order of stations in isdlite and all_vvmx_aut_60 to match the order in isd_ideam_intersect, this is, first station in column1, and last station in last column.
isdlite_stationssample = select(isdlite, paste0("X", isd_ideam_intersect["usaf" ,1:length(stationssample[,1])]))
all_vvmx_aut_60_stationssample = select(all_vvmx_aut_60, paste0("X", isd_ideam_intersect["ideam_id" ,1:length(stationssample[,1])]))


#Create one xts object for isdlite and ideam sorted data
isdlite_stationssample_xts = xts(x=as_tibble(isdlite_stationssample), order.by = timestamp_isdlite)
#summary(isdlite_stationssample_xts)
all_vvmx_aut_60_stationssample_xts = xts(x=as_tibble(all_vvmx_aut_60_stationssample), order.by = timestamp_all_vvmx_aut_60)

#ERA5 stuff
#
#
#
#
#
#
#
#
source('comparing_sources_pqrs_20199050080932_VV_AUT_ERA5_2.r')


#
#Join IDEAM and ISD to create a common plot
dev.off()
for(i in 1:length(stationssample[,1]))
  #for (i in 1:1)
{
  statisd = stationssample[i, 1]
  statisdX = paste ("X", statisd, sep="")
  #statisd = paste('"', statisd, '"', sep = "")
  #statisdX = paste('"', statisdX, '"', sep = "")

  statideam = stationssample[i, 2]
  statideamX = paste ("X", statideam, sep="")
  #statideam = paste('"', statideam, '"', sep = "")
  #statideamX = paste('"', statideamX, '"', sep = "")

  #query5 = paste("select isd_lite_unstack.mydatetime, ideam_vvmx_60.", statideam, " AS ", statideamX, ", isd_lite_unstack.",  statisd, " AS ", statisdX, " from isd_lite_unstack INNER JOIN ideam_vvmx_60 ON (isd_lite_unstack.mydatetime = ideam_vvmx_60.mydatetime) where isd_lite_unstack.mydatetime IS NOT NULL AND ", "ideam_vvmx_60.mydatetime IS NOT NULL AND ideam_vvmx_60.", statideam, " IS NOT NULL AND isd_lite_unstack.", statisd, " IS NOT NULL ORDER BY isd_lite_unstack.mydatetime", sep = "")
  currentStationIsd = select(isdlite_stationssample, statisdX)
  currentStationIsdXTS = na.omit(xts(x=as_tibble(currentStationIsd), order.by = timestamp_isdlite))

  currentStationIDEAM = select(all_vvmx_aut_60_stationssample, statideamX)
  currentStationIDEAMXTS = na.omit(xts(x=as_tibble(currentStationIDEAM), order.by = timestamp_all_vvmx_aut_60))

  isdideamera5 = merge(currentStationIDEAMXTS, currentStationIsdXTS, join='inner')
  isdideamera5 = as.data.frame(isdideamera5)
  isdideamera5$mydatetime = rownames(isdideamera5)
  isdideamera5 = select(isdideamera5, mydatetime, starts_with("X"))
  #isdideamera5 = tbl(con1, sql(query5))
  #isdideamera5 = as_tibble(select(isdideamera5, everything()))

  if (length(isdideamera5) > 0) {
    #timefilter = as_tibble(select(isdideamera5, mydatetime))$mydatetime
    my_datetime <- as.POSIXlt(as_tibble(select(isdideamera5, mydatetime))$mydatetime, format = "%Y-%m-%d %H:%M:%S")

    for (a in 1:length(isd_ideam_intersect[1,])) {
      if(isd_ideam_intersect[4,a] == stationssample[i, 2]){
        lonindex = isd_ideam_intersect[14,a]
        latindex = isd_ideam_intersect[15,a]
        xtsvar=eval(parse(text = paste("station_slice", a, isd_ideam_intersect["ideam_id",a], sep = "_")))
        xtsvar = xtsvar[my_datetime]
        era5 = as.matrix(xts(j=1, xtsvar))
        isdlon = isd_ideam_intersect["isd_lon",a]
        isdlat = isd_ideam_intersect["isd_lat",a]
        ideamlon = isd_ideam_intersect["ideam_lon",a]
        ideamlat = isd_ideam_intersect["ideam_lat",a]
        era5lonindex = isd_ideam_intersect["ideam_lon_index_era5",a]
        era5latindex = isd_ideam_intersect["ideam_lat_index_era5",a]
        era5station = isd_ideam_intersect["ideam_station_era5",a]
        colnames(era5) = c("era5")
        break
      }
    }
    isdideamera5 = cbind(isdideamera5, era5)

    #create timestamp object
    mytimestamp <- as.POSIXct(isdideamera5$mydatetime,format="%Y-%m-%d %H:%M:%S", tz="UTC")

    #Create XTS object
    assign(paste("isdideamera5", i, isd_ideam_intersect["ideam_id",i], sep = "_"), xts(x=isdideamera5[, -1], order.by = mytimestamp))

    #Plot graphics
    somePDFPath = paste(paste("../data/isdideamera5", i, statideam, sep = "_"), "pdf", sep=".")
    #pdf(file=somePDFPath,  paper="a4r", width = 0, height = 0)
    #par(mfrow = c(2,1))
    xtsvar = eval(parse(text = paste("isdideamera5", i, isd_ideam_intersect["ideam_id",i], sep = "_")))
    title = paste(
      paste("IDEAM: ", stationssample[i, 2], ". Lon: ", ideamlon, ". Lat: ", ideamlat, sep=""),
      paste("ISD: ", stationssample[i, 1], ". Lon: ", isdlon, ". Lat: ", isdlat, sep=""),
      paste("ERA5.  Col: ", era5lonindex, ". Row: ", era5latindex, ". Station: ", era5station, sep= "")
      , sep="\n")
    #par(new = FALSE)
    plot.xts(xtsvar, main = title, major.ticks="years", format.labels = "%b-%d\n%Y", add=FALSE)
    print(addLegend("top",
                    legend.names = c(paste0("IDEAM: ", stationssample[i, 2]), paste0("ISD: ", stationssample[i, 1]), "ERA5"), col=c("black", "red", "green"),
                    bg="white", bty="o", lty=c(1, 1, 1), lwd=c(2, 2, 2)))
	assign(paste0("comparison", i, "1"), recordPlot())
	saveRDS(eval(parse(text=paste0("comparison", i, "1"))), paste0("comparison", i, "1", ".rds"))

    title1 = paste(
      paste("IDEAM: ", stationssample[i, 2], ". Lon: ", ideamlon, ". Lat: ", ideamlat, sep=""),
      paste("ERA5.  Col: ", era5lonindex, ". Row: ", era5latindex, ". Station: ", era5station, sep= "")
      , sep="\n")
    par(new = FALSE)
    plot.zoo(xtsvar[,1], xtsvar[,3], xlab=paste0("IDEAM: ", stationssample[i, 2]), ylab="ERA5", main=title1)
    abline(coef=c(0,1), col="red")
	assign(paste0("comparison", i, "2"), recordPlot())
	saveRDS(eval(parse(text=paste0("comparison", i, "2"))), paste0("comparison", i, "2", ".rds"))
    title2 = paste(
      paste("ISD: ", stationssample[i, 1], ". Lon: ", isdlon, ". Lat: ", isdlat, sep=""),
      paste("ERA5.  Col: ", era5lonindex, ". Row: ", era5latindex, ". Station: ", era5station, sep= "")
      , sep="\n")
    par(new = FALSE)
    plot.zoo(xtsvar[,2], xtsvar[,3], xlab=paste0("ISD: ", stationssample[i, 1]), ylab="ERA5", main=title2)
    abline(coef=c(0,1), col="red")
	assign(paste0("comparison", i, "3"), recordPlot())
	saveRDS(eval(parse(text=paste0("comparison", i, "3"))), paste0("comparison", i, "3", ".rds"))
    title3 = paste(
      paste("IDEAM: ", stationssample[i, 2], ". Lon: ", ideamlon, ". Lat: ", ideamlat, sep=""),
      paste("ISD: ", stationssample[i, 1], ". Lon: ", isdlon, ". Lat: ", isdlat, sep="")
      , sep="\n")
    par(new = FALSE)
    plot.zoo(xtsvar[,1], xtsvar[,2], xlab=paste0("IDEAM: ", stationssample[i, 2]), ylab=paste0("ISD: ", stationssample[i, 1]),
	main=title3)
    abline(coef=c(0,1), col="red")
	assign(paste0("comparison", i, "4"), recordPlot())
	saveRDS(eval(parse(text=paste0("comparison", i, "4"))), paste0("comparison", i, "4", ".rds"))
    #dev.off()
  }
}

