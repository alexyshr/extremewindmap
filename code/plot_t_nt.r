#return level when both(t, and nt) exists
if (length(imp.vals$t.series.dt) > 0) {
  #Graphics for non-thunderstorm
  if (length(imp.vals$nt.series.dt) > 0) {
    #Intensity Function Poisson Process. NIST.SP.500-301.pdf, page 29
    Ip.t <- function(x) {1/z4 *(1+0.00001*((x-z3)/z4))^((-1/0.00001)-1)}
    Ip.nt <- function(x) {1/z7 *(1+0.00001*((x-z6)/z7))^((-1/0.00001)-1)}
    at = imp.vals$n.thunders.per.year*(1/24)
    ant = 365-at

    #Representing equation 4, page 16 (NIST.SP.500-301.pdf), but using Ip
    returnlevel_equation <- function(lower, upper, at, ant){ #lower is the value of velocity to integrate with it to infinite
      integrated.t <- integrate(Ip.t, lower=lower, upper=upper)$value
      integrated.nt <- integrate(Ip.nt, lower=lower, upper=upper)$value
      value = (at * integrated.t) + (ant * integrated.nt)
    }

    #Calculate excedence probability for velocities 1 to 600 using Poisson Process
    probabilitiesreturnlevel = sapply(yvels, returnlevel_equation, upper=Inf, at=at, ant=ant)
    veocitiesfortypicalreturnperiodsP <- approx(x=probabilitiesreturnlevel, y=yvels, xout=typicalExcedenceProbabilities)$y  #Interpolate tipical excedence probabilities to get velocities using Poisson

    zzz[zz,76:86] = veocitiesfortypicalreturnperiodsP #Velocities for typical return periods Poisson


    #
    #Plot: Page 6
    maxy = max(veocitiesfortypicalreturnperiodsP) + 0.15*max(veocitiesfortypicalreturnperiodsP)
    #Hazards Curves using Intensity Function of Poisson Process
    plot(x= 1/probabilitiesreturnlevel, y=yvels, xlab="Return Periods (Years) - Poisson Process Intensity Function",
         ylab="Velocities Km/h", main=paste("Declustered - Thunderstorms and Non-Thunderstorm - Hazard Curve - Station:", number, sep=" "),
         xlim=c(0,10000), ylim = c(0, maxy), type="l", col="blue", lwd=2)
    points(x= tipicalReturnPeriods, y= veocitiesfortypicalreturnperiodsP, col="red", pch=20)
    text(x = tipicalReturnPeriods, y = veocitiesfortypicalreturnperiodsP, labels = paste0("(",tipicalReturnPeriods,",",round(veocitiesfortypicalreturnperiodsP, digits=1),")"), cex=0.8, pos = 4)
    mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
    legend("bottomright", c("Hazard Curve"),
           bty = "n",
           col = "blue",
           lty = c(1), lwd = c(2),
           pch = c(NA),
           pt.bg = c(NA),
           pt.cex = c(0))
    box(lty = 1, col = 'black', lwd=0.5)
    assign(paste0("myprint", numberofplots), recordPlot())
    #saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0("myprint", numberofplots, ".rds"))
    numberofplots = numberofplots + 1
  }
}
