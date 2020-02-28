#Create Raw Data Statistics and Send to CSV - Non Thunderstorm
if (length(raw.data.nt$date.time) > 0) {
  years = generate_stats_time_serie(raw.data.nt, "speed.kph", raw.data.nt$date.time, "years")
  months = generate_stats_time_serie(raw.data.nt, "speed.kph", raw.data.nt$date.time, "months")
  weeks = generate_stats_time_serie(raw.data.nt, "speed.kph", raw.data.nt$date.time, "weeks")

  write.xlsx(years, file=statsfile, sheetName="nt_years", append=TRUE, row.names=TRUE)
  write.xlsx(months, file=statsfile, sheetName="nt_months", append=TRUE, row.names=TRUE)
  write.xlsx(weeks, file=statsfile, sheetName="nt_weeks", append=TRUE, row.names=TRUE)

  #Search time differences in days between consecutive samples greather than threshold in days (last parameter next function)
  thresholdindays = 30
  holesindays = locate_holes_time_serie(raw.data.nt, "speed.kph", raw.data.nt$date.time, thresholdindays)
  if (length(holesindays) == 0){
    holesindays = "No holes!"
  }
  write.xlsx(holesindays, file=statsfile, sheetName=paste0("nt_gaps",thresholdindays,"days"), append=TRUE, row.names=FALSE)
  #Plot time serie
  opar <- par(no.readonly = TRUE)
  print(plotxts(data=raw.data.nt, variable="speed.kph", time=raw.data.nt$date.time, cex.main=0.2, major.ticks="years",
                xlab=paste0("Page ",numberofplots, " - Time Series Plot for Non-Thunderstorm ('nt') - Station: ", number),
                main = paste0("Station ID: ",  number, "\nWind Velocity [Km/h]")))
  assign(paste0("myprint", numberofplots), recordPlot())
  saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0(outputpath, "myprint", numberofplots, ".rds"))
  numberofplots = numberofplots + 1
  #dev.off()
  par(opar)
}
