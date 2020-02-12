if (length(raw.data.tibble$date.time) > 0) {

  years = generate_stats_time_serie(raw.data.tibble, "speed.kph", raw.data.tibble$date.time, "years")
  months = generate_stats_time_serie(raw.data.tibble, "speed.kph", raw.data.tibble$date.time, "months")
  weeks = generate_stats_time_serie(raw.data.tibble, "speed.kph", raw.data.tibble$date.time, "weeks")

  #myxts = na.omit(xts(x=select(raw.data.tibble, "speed.kph"), order.by = raw.data.tibble$date.time))
  #Split dataset by year
  #xts5_yearly <- split(myxts,f="years") #Convert to [[]]
  #sc = lapply(xts5_yearly, cumsum) #Calcular suma cumulativa but still [[]]
  #do.call(rbind,  sc) #Unstack [[]] and leave all in one index []

  #statistics <- list()
  #statistics[[1]] = years
  #statistics[[2]] = months
  #statistics[[3]] = weeks
  #lapply(statistics, write, statsfile, append=TRUE)

  library(xlsx)
  write.xlsx(years, file=statsfile, sheetName="all_years", row.names=TRUE)
  write.xlsx(months, file=statsfile, sheetName="all_months", append=TRUE, row.names=TRUE)
  write.xlsx(weeks, file=statsfile, sheetName="all_weeks", append=TRUE, row.names=TRUE)

  #Search time differences in days between consecutive samples greather than threshold in days (last parameter next function)
  thresholdindays = 30
  holesindays = locate_holes_time_serie(raw.data.tibble, "speed.kph", raw.data.tibble$date.time, thresholdindays)
  if (length(holesindays) == 0){
    holesindays = "No holes!"
  }
  write.xlsx(holesindays, file=statsfile, sheetName=paste0("all_gaps",thresholdindays,"days"), append=TRUE, row.names=FALSE)  

  #Plot time serie
  #library(xts)
  #myxts = na.omit(xts(x=select(raw.data.tibble, "speed.kph"), order.by = raw.data.tibble$date.time))
  #par(oma = c(2,0,0,0))
  #main=paste0("Time Series Plot for Raw.Data\nStation: ", number, " - Wind Velocity [Km/h]")
  #print(plot.xts(myxts, main=main, major.ticks="year", format.labels = "%b-%d\n%Y",
  #         col="green", legend.loc = "top", cex.main=0.2))
  #mtext(side = 1, text = paste0("Page ",numberofplots, " - Time Series Plot for Raw.Data - Station: ", number), outer = TRUE)
  print(plotxts(data=raw.data.tibble, variable="speed.kph", time=raw.data.tibble$date.time,
                cex.main=0.2, major.ticks="years",
                xlab=paste0("Page ",numberofplots, " - Time Series Plot for Raw.Data - Station: ", number),
                main = paste0("Station ID: ",  number, "\nWind Velocity [Km/h]")))
  assign(paste0("myprint", numberofplots), recordPlot())
  saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0(outputpath, "myprint", numberofplots, ".rds"))				
  numberofplots = numberofplots + 1
}
