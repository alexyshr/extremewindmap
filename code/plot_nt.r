#Graphics for non-thunderstorm
if (length(imp.vals$nt.series.dt) > 0) {
  #___
  #___a) Poisson Process using Shape parameter!!
  #___

  #Intensity Function Poisson Process. NIST.SP.500-301.pdf, page 29
  Ip.nt <- function(x) {1/z7 *(1+0.00000001*((x-z6)/z7))^((-1/0.00000001)-1)}
  #Ip.nt <- function(x) {1/z7 * exp(-(x-z6)/z7)}

  ant = 365 #typical amount of non-thunderstorm in a year. As "at" = 0, then ant = 365 -0
  if (length(imp.vals$t.series.dt) > 0) {
    ant = 365 - (imp.vals$n.thunders.per.year*(1/24)) #But if nt is not zero i have to substract at it
  }

  #Representing equation 4, page 16 (NIST.SP.500-301.pdf), but using Ip
  returnlevel_equation <- function(lower, upper, ant){ #lower is the value of velocity to integrate with it to infinite
    integrated.nt <- integrate(Ip.nt, lower=lower, upper=upper)$value
    value = ant * integrated.nt
  }

  #Calculate excedence probability for velocities 1 to 600 using Poisson Process
  paP = sapply(yvels, returnlevel_equation, upper=Inf, ant=ant)

  #Probability density according to NIST.SP.500-301.pdf. Equation after equation 2, page 15
  densityPoissonProcess <- function (velocity, series, threshold){#, time){  #lower is the value of velocity to integrate with it up to infinite
    y1 = min(series)
    y2 = max(series)
    #  t1 = min(serie.dt)
    #  t2 = max(serie.dt)
    #t2_minus_t1 <- as.numeric(difftime(time[length(time)],
    #                        time[1],
    #                       units="days"))
    integrand <- Ip.nt
    #integrated <- integrate(integrand, lower=y1, upper=y2)$value
    integrated <- integrate(integrand, lower=threshold, upper=y2)$value

    #integrated <- integrate(integrand, lower=lower, upper=upper)$value
    #density <- Ip.nt(velocity)/(t2_minus_t1*integrated)  # Equation after equation 2, page 15
    density <- Ip.nt(velocity)/(integrated)  # Equation after equation 2, page 15
  }

  yv = seq(from=min(imp.vals$nt.series), to= max(imp.vals$nt.series), length.out=1000)
  dfP = sapply(yv, densityPoissonProcess, series=imp.vals$nt.series, threshold=z5)#, time=imp.vals$nt.series.dt)


  #Plot: Page 2: Density Function from Intensity Function of Poisson Process
  #Plot dfP - probability density function using intensity function of Poisson Process
  library(RColorBrewer)
  cols <- brewer.pal(9,"Set1")

  plot(x= yv, y=dfP, xlab="Declustered - Non-Thunderstorm - Velocities [Km/h]", ylab="Density Function - Poisson Process",
       main=paste("Density Function from Intensity Function of Poisson Process\nStation:", number, sep=" "),
       type="l", col=cols[3], lwd=4)
  mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  hist(imp.vals$nt.series, probability = TRUE, add=TRUE, col="cadetblue3")
  lines(x= yv, y=dfP, col=cols[3], lwd=4)
  curve(dgpd_pp(x, location=z6, shape=0.00000001, scale=z7, threshold=z5), col=cols[5], add=TRUE, lwd=2, lty= 2)

  legend("topright", c("Data-Empirical", "Fitted-Theoretical POT-Poisson Process", "Fitted-Theoretical POT-GPD Equivalence"),
         bty = "n",
         col = c("cadetblue3", cols[3], cols[5]),
         lty = c(0,1,2), lwd = c(0,4,2),
         pch = c(22, NA, NA),
         pt.bg = c("cadetblue3", NA, NA),
         pt.cex = 2)
  assign(paste0("myprint", numberofplots), recordPlot())
  #saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0("myprint", numberofplots, ".rds"))
  numberofplots = numberofplots + 1


  #Fit Pareto of Poisson Process to Get Graphics and Goodness of Fit
  gpdScale = (-0.00001)*((-z7/0.00001)+z6-z5)
  gpdLocation = get("z5")
  attr(gpdLocation, "names") = NULL
  param5 = list(scale = gpdScale, location=gpdLocation, shape=0.00001)
  param6 = list(dummy = 1)
  #Fit equivalent GPD
  fpnt <- fitdistrplus::fitdist(imp.vals$nt.series, "mygpd",
                                start=param6, fix.arg=param5, lower = c(0,0),
                                method="mle", optim.method = "L-BFGS-B")
  #print(fp)
  #summary(fp)
  #cdfcomp(fp)
  #denscomp(fp)
  #ppcomp(fp)
  #qqcomp(fp)
  #quantile(fp, probs = c(0.05, 0.1, 0.2))
  mygof= fitdistrplus::gofstat(fpnt, discrete=FALSE) # goodness-of-fit statistics
  mykstest = stats::ks.test(imp.vals$nt.series, y="pmygpd", location=gpdLocation, scale=gpdScale, shape=0.00001, dummy=1)

  par(cex=1.2, bg="white")
  plot(fpnt, lwd=2, col="cadetblue3")

  mtext(side = 1, text = paste0("Page ", numberofplots, " - Declustered - Non-Thunderstorm - POT-GPD Equivalent. Location: ",
                                round(fpnt$fix.arg$location, digits = 2), ". Scale: ", round(fpnt$fix.arg$scale, digits = 2),
                                ". Shape: ", round(fpnt$fix.arg$shape, digits = 2)), outer = TRUE)
  assign(paste0("myprint", numberofplots), recordPlot())
  #saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0("myprint", numberofplots, ".rds"))
  numberofplots = numberofplots + 1

  #RMSE
  mx <- seq(min(imp.vals$nt.series), max(imp.vals$nt.series), length = 10000)
  empCDF <- stats::ecdf(imp.vals$nt.series)
  #D of Kolmogorov Smirnov
  #max(abs(empCDF(mx) - pmygpd(mx, scale = fpnt$fix.arg$scale, location =fpnt$fix.arg$location, shape = 0.00001, dummy=1)))
  #rmse
  rmse = sqrt(mean((empCDF(mx) - pmygpd(mx, scale = fpnt$fix.arg$scale, location =fpnt$fix.arg$location, shape = 0.00001, dummy=1))^2))

  myfitgpd = list("Fit_fitdistrplus::fitdist" = capture.output(fpnt),
                  "GoodnessOfFit_fitdistrplus::gofstat"=capture.output(mygof),
                  "KolmogorovSmirnovTest_stats::ks.test"=capture.output(mykstest),
                  "RMSE"=capture.output(rmse))
  write.xlsx(capture.output(myfitgpd), file=fnfitted, sheetName="nt_POT-GPD-Equivalent", append=TRUE, row.names=TRUE)


  #Cumulative distribution function according to NIST.SP.500-301.pdf. Equation 5, page 17
  distributionPoissonProcess <- function(velocity, threshold){  #lower is the value of velocity to integrate with it up to infinite
    integrand <- Ip.nt
    integrated_tresh_to_vel <- integrate(integrand, lower=threshold, upper=velocity)$value
    integrated_tresh_to_inf <- integrate(integrand, lower=threshold, upper=Inf)$value
    #distribution <- 1- (integrated_tresh_to_vel/integrated_tresh_to_inf)#NIST.SP.500-301.pdf. Equation 5, page 17
    #Be ware I am not using 1 - () as is stated in equation 5
    distribution <- (integrated_tresh_to_vel/integrated_tresh_to_inf)#NIST.SP.500-301.pdf. Equation 5, page 17
  }

  #Calculate cumulative distribution function values for velocities y1 to y2 (Domain PP) using Poisson Process
  DfP = sapply(yv, distributionPoissonProcess, threshold=z5)


  #Plot: Page 3 - Cumulative distribution function with intensity function of Poisson Process
  #Plot Dfp- Cumulative distribution function  with intensity funtion of Poisson Process
  library(RColorBrewer)
  cols <- brewer.pal(9,"Set1")
  plot(x= yv, y=DfP, xlab="Declustered - Non-Thuderstorm - Velocities [Km/h]", ylab="Distribution Function - Poisson Process",
       main=paste("Cumulative Distribution Function from Intensity Function of Poisson Process\nStation:", number, sep=" "),
       type="l", col=cols[3], lwd=4)
  mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  plot(stats::ecdf(imp.vals$nt.series), col="cadetblue3", add=TRUE)
  lines(x= yv, y=DfP, col=cols[3], lwd=4)
  curve(pgpd_pp(x, location=z6, shape=0.00000001, scale=z7, threshold=z5), col=cols[5], add=TRUE, lwd=2, lty= 2)
  legend("bottomright", c("Data-Empirical", "Fitted-Theoretical POT-Poisson Process", "Fitted-Theoretical POT-GPD Equivalence"),
         bty = "o",
         col = c("cadetblue3", cols[3], cols[5]),
         lty = c(1, 1, 2), lwd = c(1, 4, 2),
         pch = c(21, NA, NA),
         pt.bg = c("cadetblue3", NA, NA),
         pt.cex = 1, box.lwd = 0, box.col = "white", bg = "white")
  box(lty = 1, col = 'black', lwd=0.5)
  box(lty = 1, col = 'black', lwd=0.5)
  assign(paste0("myprint", numberofplots), recordPlot())
  #saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0("myprint", numberofplots, ".rds"))
  numberofplots = numberofplots + 1

  #Plot: Page 4
  #Histogram and fitted density function from density function of NIST.SP.500-301.pdf. Equation after equation 2, page 15
  #Whitout Substracting Threshold
  overthreshold = imp.vals$nt.series >= z5 #Be sure only data above threshold are included
  #hist(imp.vals$nt.series[overthreshold], prob=TRUE, main=paste("Fitted density function from Intensity Function of Poisson Process\nWhitout Substracting Threshold\nStation:", number, sep=" "))
  #lines(x=yvels, y=dfP, col="red")
  #mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  #numberofplots = numberofplots + 1

  #Plot: Page 5
  #Histogram and fitted density function from density function of NIST.SP.500-301.pdf. Equation after equation 2, page 15
  #Subtracting Threshold
  #overthreshold = imp.vals$nt.series >= z5 #Be sure only data above threshold are included
  #hist(imp.vals$nt.series[overthreshold]-nt.thresh, prob=TRUE, main=paste("Fitted density function from Intensity Function of Poisson Process\nSubstracting Threshold\nStation:", number, sep=" "))
  #lines(x=yvels, y=dfP, col="red")
  #mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  #numberofplots = numberofplots + 1

  #Calculate velocities for typical return periods using Poisson Process
  veocitiesfortypicalreturnperiodsP <- approx(x=paP, y=yvels, xout = typicalExcedenceProbabilities)$y  #Interpolate tipical excedence probabilities to get velocities using Poisson

  zzz[zz,43:53] = veocitiesfortypicalreturnperiodsP #Velocities for typical return periods Poisson

  #tipicalReturnPeriods = c("10","20","50","100","250","500","700","1000","1700","3000","7000")
  #
  #Plot: Page 6
  #Hazards Curves using Intensity Function of Poisson Process
  library(RColorBrewer)
  cols <- brewer.pal(9,"Set1")
  maxy = max(veocitiesfortypicalreturnperiodsP) + 0.15*max(veocitiesfortypicalreturnperiodsP)
  plot(x= 1/paP, y=yvels, xlab="Return Periods (Years) - Poisson Process Intensity Function",
       ylab="Velocities Km/h", main=paste("Declustered - Non-Thunderstorms - Hazard Curve - Station:", number, sep=" "),
       xlim=c(0,10000), ylim=c(0, maxy), type="l", col=cols[3], lwd=3)
  myx= 1:10000
  myrlpp = sapply(myx, returnlevelPP, location=z6, shape=0.00001, scale=z7)
  lines(myx, myrlpp, col=cols[4], lwd=2, lty=5)

  #curve (returnlevelGPD(x, location=z3, shape=0.00001, scale=z4, threshold=z2), col="black", add=TRUE, lwd=2)
  myrlgpd2 = sapply(myx, returnlevelGPD, location=z6, shape=0.00001, scale=z7, threshold=z5)
  lines(myx, myrlgpd2, col=cols[5], lwd=1, lty=2)

  points(x= tipicalReturnPeriods, y= veocitiesfortypicalreturnperiodsP, col=cols[9], pch=20)
  text(x = tipicalReturnPeriods, y = veocitiesfortypicalreturnperiodsP, labels = paste0("(",tipicalReturnPeriods,",",round(veocitiesfortypicalreturnperiodsP, digits=1),")"), cex=0.8, pos = 4)
  mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  legend("bottomright", c("POT-PP (Intensity Function Integral)","POT-PP (formula)", "POT-GPD equivalence"),
         bty = "n",
         col = cols[3:5],
         lty = c(1, 5, 2), lwd = c(3, 2, 1),
         pch = c(NA, NA, NA),
         pt.bg = c(NA, NA, NA),
         pt.cex = c(0,0,0))
  box(lty = 1, col = 'black', lwd=0.5)
  assign(paste0("myprint", numberofplots), recordPlot())
  #saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0("myprint", numberofplots, ".rds"))
  numberofplots = numberofplots + 1


  #___
  #___b) Poisson Process - Gumbel like tile length. Shape = 0
  #___

  #Intensity Function Poisson Process - Gumbel like tile length. Shape = 0. NIST.SP.500-301.pdf, page 32, equation 10
  Ig.nt <- function(x) {1/z7 * exp(-(x-z6)/z7)}
  ant = 365 #typical amount of non-thunderstorm in a year. As "at" = 0, then ant = 365 -0
  if (length(imp.vals$t.series.dt) > 0) {
    ant = 365 - (imp.vals$n.thunders.per.year*(1/24)) #But if nt is not zero
  }

  #Representing equation 4, page 16 (NIST.SP.500-301.pdf), but using Ip
  returnlevel_equation_gumbel <- function(lower, upper, ant){  #lower is the value of velocity to integrate with it up to infinite
    integrated.nt <- integrate(Ig.nt, lower=lower, upper=upper)$value
    value = ant * integrated.nt
  }

  #Calculate excedence probability for velocities 1 to 600 using Poisson Process - Gumbel Intensity Function
  paG = sapply(yvels, returnlevel_equation_gumbel, upper=Inf, ant=ant)

  #Calculate velocities for typical return periods using Poisson Process - Gumbel Intensity Function
  veocitiesfortypicalreturnperiodsG <- approx(x=paG, y=yvels, xout = typicalExcedenceProbabilities)$y  #Interpolate tipical excedence probabilities to get velocities using Gumbel
  zzz[zz,54:64] = veocitiesfortypicalreturnperiodsG #Velocities for typical return periods Gumbel Tile Length Intensity Function

  #Plot: Page 7
  maxy = max(veocitiesfortypicalreturnperiodsG) + 0.15*max(veocitiesfortypicalreturnperiodsG)

  library(RColorBrewer)
  cols <- brewer.pal(9,"Set1")

  #Hazards Curves using Intensity Function of Gumbel like tail length
  plot(x= 1/paG, y=yvels, xlab="Return Periods (Years) - Gumbel like tail Intensity Function of Poisson Process",
       ylab="Velocities Km/h", main=paste("Declustered - Non-Thunderstorms - Hazard Curve - Station:", number, sep=" "),
       xlim=c(0,10000), ylim=c(0,maxy), type="l", col=cols[3], lwd=3)
  myx= 1:10000
  myrlpp = sapply(myx, returnlevelPP, location=z6, shape=0.00001, scale=z7)
  lines(myx, myrlpp, col=cols[4], lwd=2, lty=5)

  #curve (returnlevelGPD(x, location=z6, shape=0.00001, scale=z7, threshold=z7), col="black", add=TRUE, lwd=2)
  myrlgpd2 = sapply(myx, returnlevelGPD, location=z6, shape=0.00001, scale=z7, threshold=z5)
  lines(myx, myrlgpd2, col=cols[5], lwd=1, lty=2)

  points(x= tipicalReturnPeriods, y= veocitiesfortypicalreturnperiodsG, col=cols[9], pch=20)
  text(x = tipicalReturnPeriods, y = veocitiesfortypicalreturnperiodsG, labels = paste0("(",tipicalReturnPeriods,",",round(veocitiesfortypicalreturnperiodsG, digits=1),")"), cex=0.8, pos = 4)
  mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  legend("bottomright", c("POT-PP (Gumbel like - Intensity Function Integral)","POT-PP (formula)", "POT-GPD equivalence"),
         bty = "n",
         col = cols[3:5],
         lty = c(1, 5, 2), lwd = c(3, 2, 1),
         pch = c(NA, NA, NA),
         pt.bg = c(NA, NA, NA),
         pt.cex = c(0,0,0))
  box(lty = 1, col = 'black', lwd=0.5)
  assign(paste0("myprint", numberofplots), recordPlot())
  #saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0("myprint", numberofplots, ".rds"))

  numberofplots = numberofplots + 1



  #___
  #___c) Using Gumbel Distribution replacing fitted parameters (location=z6, scale=z7)
  #___



  #Probability density function using Gumbel distribution
  #Plot: Page 8

  library(RcmdrMisc)
  .x <- seq(-300, 200, length.out=1000)
  dfG = RcmdrMisc::dgumbel(.x, location=z6, scale=z7)
  plot(.x, dfG, col="red",
       xlab="Non-Thunderstorms - Velocities Km/h", ylab="Density Function - Gumbel Distribution",
       main=paste("Gumbel Density Function, but using parameters of Poisson Process\n", "Location=", round(z6,2), " Scale=", round(z7,2), "\n", "Station:", number),
       type="l", lwd=1, ylim=c(0,0.08))
  #plotDistr(.x, dgumbel(.x, location=0, scale=1), cdf=FALSE, xlab="x",
  #          ylab="Density", main=paste("Gumbel Distribution:  Location=0, Scale=1"))
  mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  hist(imp.vals$nt.series[overthreshold], prob=TRUE, add=T, col="cadetblue3")
  pdfG <- function(x) {1/z7 *exp(-(x-z6)/z7)*exp(-exp(-(x-z6)/z7))}
  #curve(pdfG, add=T)
  lines(.x, dfG, col="red")
  legend("topright", c("Data-Empirical", "Fitted-Theoretical"),
         bty = "n",
         col = c("cadetblue3", "red"),
         lty = c(0, 1), lwd = c(0, 1),
         pch = c(22, NA),
         pt.bg = c("cadetblue3", NA),
         pt.cex = 2)
  assign(paste0("myprint", numberofplots), recordPlot())
  #saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0("myprint", numberofplots, ".rds"))

  numberofplots = numberofplots + 1

  #Draw both pdf. Alexys H
  #pdfGfgev <- function(x) {1/scalebyfgev *exp(-(x-locbyfgev)/scalebyfgev)*exp(-exp(-(x-locbyfgev)/scalebyfgev))}
  #curve(pdfGfgev, add=T, col="cadetblue3")


  #pdfG <- function(x) {1/z7 *exp(-(x-0)/z7)*exp(-exp(-(x-0)/z7))}
  #curve(pdfG, add=T)


  #Cumulative distribution function with Gumbel distribution
  #Plot: Page 9
  .x <- seq(-300, 200, length.out=1000)
  DfG = RcmdrMisc::pgumbel(.x, loc=z6, scale=z7)
  plot(.x, DfG, lwd=2,
       xlab="Non-Thunderstorms - Velocities Km/h", ylab="Cumulative Distribution Function - Gumbel Distribution",
       main=paste("Gumbel Cumulative Distribution, but using parameters of Poisson Process \n", "Location=", round(z6,2), " Scale=", round(z7,2), "\n", "Station:", number),
       type="l", col="red")
  #lines(yvels, DfG, col="red", lty=4)
  mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  plot(stats::ecdf(imp.vals$nt.series), col="cadetblue3", add=TRUE)
  lines(.x, DfG, col="red")
  legend("bottomright", c("Data-Empirical", "Fitted-Theoretical"),
         bty = "n",
         col = c("cadetblue3", "red"),
         lty = c(1, 1), lwd = c(1, 2),
         pch = c(21, NA),
         pt.bg = c("cadetblue3", NA),
         pt.cex = 1)
  assign(paste0("myprint", numberofplots), recordPlot())
  #saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0("myprint", numberofplots, ".rds"))

  numberofplots = numberofplots + 1

  #Plot: Page 10
  #Histogram and fitted density function from Gumbel distribution
  overthreshold = imp.vals$nt.series >= z5 #Be sure only data above threshold are included
  hist(imp.vals$nt.series[overthreshold], xlab="Declustered Non-Thunderstorm", prob=TRUE, col="cadetblue3",
       main=paste("Fitted Gumbel density function using parameters of Poisson Process\n", "Location=", round(z6,2), " Scale=", round(z7,2), "\n", "Station:", number))
  curve(RcmdrMisc::dgumbel(x, location=z6, scale=z7), add=TRUE, col="red", lwd=2)
  #lines(x=yvels, y=dfP, col="green")
  #lines(densityPoissonProcess(lower=x, upper=inf), add = TRUE)
  mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  legend("topright", c("Data-Empirical", "Fitted-Theoretical"),
         bty = "n",
         col = c("cadetblue3", "red"),
         lty = c(0, 1), lwd = c(0, 2),
         pch = c(22, NA),
         pt.bg = c("cadetblue3", NA),
         pt.cex = 2)
  assign(paste0("myprint", numberofplots), recordPlot())
  #saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0("myprint", numberofplots, ".rds"))

  numberofplots = numberofplots + 1

  #Following section "Using the Fitted Model" of Extreme Wind Speeds: Overview (https://www.itl.nist.gov/div898/winds/overview.htm)
  #Calculation of Velocities for Return Periods using Percent Point Function of Gumbel (Quantiles: qgumbel)

  gumbelVelocitiesQuantileFunction <- function(mri){  #mri: return intervals from 1 to 3000
    RcmdrMisc::qgumbel((1-(1/mri)), location=z6, scale=z7)
  }

  gumbelVelocitiesQuantileFunctionLambda <- function(mri, npy, numberofsamples, numberofsamplesoverthreshold){  #mri: return intervals from 1 to 3000
    #npy #Number of observations per year
    #w = length(imp.vals$nt.series)/npy  #Fitting period: Total observations divided in npy
    w = numberofsamples/npy
    #Observations over threshold
    #overthresh = imp.vals$nt.series > nt.thresh
    #overthreshold = imp.vals$t.series >= z2
    #lambda = length(imp.vals$nt.series[overthresh])/w
    lambda = numberofsamplesoverthreshold/w
    RcmdrMisc::qgumbel((1-(1/(lambda*mri))), location=z6, scale=z7)
  }

  gumbelVelocitiesQuantileFunctionNpy <- function(mri, npy){  #mri: return intervals from 1 to 3000
    RcmdrMisc::qgumbel((1-(1/(npy*mri))), location=z6, scale=z7)
  }

  gumbelVelocitiesQuantileFunctionBerry <- function(mri, npy, truncate){  #mri: return intervals from 1 to 3000
    probs = 1 -(1/(mri*npy))
    probs2 = (probs-truncate)/(1-truncate)
    RcmdrMisc::qgumbel(probs2, location=z6, scale=z7)
  }


  #Calculation of Velocities for Typical Return Periods using Quantile Function of Gumbel Distribution
  #veocitiesfortypicalreturnperiodsQG = sapply(tipicalReturnPeriods, function(x){qgumbel((1-(1/x)), location=z6, scale=z7)})
  veocitiesfortypicalreturnperiodsQG = sapply(tipicalReturnPeriods, gumbelVelocitiesQuantileFunction)

  #Next option using Lambda
  #npy=imp.vals$n.nthunders.per.year
  #veocitiesfortypicalreturnperiodsQG = sapply(tipicalReturnPeriods, gumbelVelocitiesQuantileFunctionLambda,
  #                                            npy=npy, numberofsamples= length(imp.vals$nt.series),
  #                                            numberofsamplesoverthreshold= length(imp.vals$nt.series[imp.vals$nt.series > nt.thresh]) )

  #This option using npy
  #npy=imp.vals$n.nthunders.per.year
  #veocitiesfortypicalreturnperiodsQG = sapply(tipicalReturnPeriods, gumbelVelocitiesQuantileFunctionNpy, npy=npy)


  #This option using berry
  #npy=imp.vals$n.nthunders.per.year
  #truncate = 1 - (sum(overthreshold)/length(imp.vals$nt.series))
  #veocitiesfortypicalreturnperiodsQG = sapply(tipicalReturnPeriods, gumbelVelocitiesQuantileFunctionBerry, npy=npy, truncate=truncate)


  zzz[zz,65:75] = veocitiesfortypicalreturnperiodsQG #Velocities for typical return periods using Gumbel quantile function

  #Calculation of velocities for return periods from 1 to 3000
  mri = 1:10000
  allVelGumbelQuantileFunction = sapply(mri, gumbelVelocitiesQuantileFunction)
  library(RColorBrewer)
  cols <- brewer.pal(9,"Set1")
  #Hazard curve using quantile function of Gumbel distribution
  #Plot: Page 11
  plot(x= mri, y=allVelGumbelQuantileFunction, xlab="Return Periods (Years) - Gumbel Quantile Function using parameters of Poisson Process",
       ylab="Velocities Km/h", main=paste("Declustered - Non-Thunderstorms - Hazard Curve - Station:", number, sep=" "),
       xlim=c(0,10000), type="l", col=cols[3], lwd=3)
  myx= 1:10000
  myrlpp = sapply(myx, returnlevelPP, location=z6, shape=0.00001, scale=z7)
  lines(myx, myrlpp, col=cols[4], lwd=2, lty=5)

  #curve (returnlevelGPD(x, location=z6, shape=0.00001, scale=z7, threshold=z7), col="black", add=TRUE, lwd=2)
  myrlgpd2 = sapply(myx, returnlevelGPD, location=z6, shape=0.00001, scale=z7, threshold=z5)
  lines(myx, myrlgpd2, col=cols[5], lwd=1, lty=2)

  points(x= tipicalReturnPeriods, y= veocitiesfortypicalreturnperiodsQG, col=cols[9], pch=20)
  text(x = tipicalReturnPeriods, y = veocitiesfortypicalreturnperiodsQG, labels = paste0("(",tipicalReturnPeriods,",",round(veocitiesfortypicalreturnperiodsQG, digits=1),")"), cex=0.8, pos = 4)
  legend("bottomright", c("Gumbel Quantile","POT-PP (formula)", "POT-GPD equivalence"),
         bty = "n",
         col = cols[3:5],
         lty = c(1, 5, 2), lwd = c(3, 2, 1),
         pch = c(NA, NA, NA),
         pt.bg = c(NA, NA, NA),
         pt.cex = c(0,0,0))
  box(lty = 1, col = 'black', lwd=0.5)
  mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  assign(paste0("myprint", numberofplots), recordPlot())
  #saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0("myprint", numberofplots, ".rds"))

  numberofplots = numberofplots + 1
}
