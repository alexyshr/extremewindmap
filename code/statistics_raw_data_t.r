#Create Raw Data Statistics and Send to CSV - Thunderstorm
if (length(raw.data.t$date.time) > 0) {
  years = generate_stats_time_serie(raw.data.t, "speed.kph", raw.data.t$date.time, "years")
  months = generate_stats_time_serie(raw.data.t, "speed.kph", raw.data.t$date.time, "months")
  weeks = generate_stats_time_serie(raw.data.t, "speed.kph", raw.data.t$date.time, "weeks")

  write.xlsx(years, file=statsfile, sheetName="t_years", append=TRUE, row.names=TRUE)
  write.xlsx(months, file=statsfile, sheetName="t_months", append=TRUE, row.names=TRUE)
  write.xlsx(weeks, file=statsfile, sheetName="t_weeks", append=TRUE, row.names=TRUE)

  #Search time differences in days between consecutive samples greather than threshold in days (last parameter next function)
  thresholdindays = 30
  holesindays = locate_holes_time_serie(raw.data.t, "speed.kph", raw.data.t$date.time, thresholdindays)

  write.xlsx(holesindays, file=statsfile, sheetName=paste0("t_gaps",thresholdindays,"days"), append=TRUE, row.names=FALSE)
  #Plot time serie
  print(plotxts(data=raw.data.t, variable="speed.kph", time=raw.data.t$date.time, cex.main=0.2, major.ticks="years",
                xlab=paste0("Page ",numberofplots," - Time Series Plot for Thunderstorm ('t') - Station: ", number),
                main = paste0("Station ID: ",  number, "\nWind Velocity [Km/h]")))
  assign(paste0("myprint", numberofplots), recordPlot())
  #saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0("myprint", numberofplots, ".rds"))
  numberofplots = numberofplots + 1
}
