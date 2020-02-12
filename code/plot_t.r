
#Graphics for thunderstorm
if (length(imp.vals$t.series.dt) > 0) {
  #___
  #___a) Poisson Process using Shape parameter!!
  #___
  
  
  #Intensity Function Poisson Process. NIST.SP.500-301.pdf, page 29
  Ip.t <- function(x) {1/z4 *(1+0.00001*((x-z3)/z4))^((-1/0.00001)-1)}
  #ans = 365-as  #typical amount of non-thunderstorm in a year. As "as" = 0, then ans = 365 -0
  at = imp.vals$n.thunders.per.year*(1/24)
  
  #Representing equation 4, page 16 (NIST.SP.500-301.pdf), but using Ip
  returnlevel_equation <- function(lower, upper, at){ #lower is the value of velocity to integrate with it to infinite
    integrated.t <- integrate(Ip.t, lower=lower, upper=upper)$value
    value = at * integrated.t
  }
  
  #Calculate excedence probability for velocities 1 to 600 using Poisson Process
  paP = sapply(yvels, returnlevel_equation, upper=Inf, at=at)
  
  #Probability density according to NIST.SP.500-301.pdf. Equation after equation 2, page 15
  densityPoissonProcess <- function (velocity, series, threshold){  #lower is the value of velocity to integrate with it up to infinite
    y1 = min(series)
    y2 = max(series)
    #  t1 = min(serie.dt)
    #  t2 = max(serie.dt)
    integrand <- Ip.t
    #integrated <- integrate(integrand, lower=y1, upper=y2)$value
    integrated <- integrate(integrand, lower=threshold, upper=y2)$value
    #integrated <- integrate(integrand, lower=lower, upper=upper)$value
    density <- Ip.t(velocity)/integrated  # Equation after equation 2, page 15
  }
  
  #Calculate density function values for velocities y1 to y1 (Domain in Velocity) using Poisson Process
  #dfP = sapply(yvels, densityPoissonProcess, upper=Inf)
  yv = seq(from=min(imp.vals$t.series), to= max(imp.vals$t.series), length.out=1000)
  dfP = sapply(yv, densityPoissonProcess, series=imp.vals$t.series, threshold=z2)
  library(RColorBrewer)
  cols <- brewer.pal(9,"Set1")
  #Plot: Page 2: Density Function from Intensity Function of Poisson Process
  #Plot dfP - probability density function using intensity function of Poisson Process
  plot(x= yv, y=dfP, xlab="Declustered - Thunderstorm - Velocities [Km/h]", ylab="Density Function - Poisson Process",
       main=paste("Density Function from Intensity Function of Poisson Process\n", "Location=", round(z3,2), " Scale=", round(z4,2), "\n", "Station:", number),
       type="l", col=cols[3], lwd=4)
  mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  hist(imp.vals$t.series, probability = TRUE, add=TRUE, col="cadetblue3")
  lines(x= yv, y=dfP, col=cols[3], lwd=4)
  curve(dgpd_pp(x, location=z3, shape=0.00000001, scale=z4, threshold=z2), col=cols[5], add=TRUE, lwd=2, lty= 2)
  legend("topright", c("Data-Empirical", "Fitted-Theoretical POT-Poisson Process", "Fitted-Theoretical POT-GPD Equivalence"),
         bty = "n",
         col = c("cadetblue3", cols[3], cols[5]),
         lty = c(0,1,2), lwd = c(0,4,2),
         pch = c(22, NA, NA),
         pt.bg = c("cadetblue3", NA, NA),
         pt.cex = 2)
  assign(paste0("myprint", numberofplots), recordPlot())
  saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0(outputpath, "myprint", numberofplots, ".rds"))		 
  numberofplots = numberofplots + 1
  
  #Fit Pareto of Poisson Process to Get Graphics and Goodness of Fit
  gpdScale = (-0.00001)*((-z4/0.00001)+z3-z2)
  gpdLocation = get("z2")
  attr(gpdLocation, "names") = NULL
  param5 = list(scale = gpdScale, location=gpdLocation, shape=0.00001)
  param6 = list(dummy = 1)
  #Fit equivalent GPD
  fp <- fitdistrplus::fitdist(imp.vals$t.series, "mygpd",
                              start=param6, fix.arg=param5, lower = c(0,0),
                              method="mle", optim.method = "L-BFGS-B")
  #print(fp)
  #summary(fp)
  #cdfcomp(fp)
  #denscomp(fp)
  #ppcomp(fp)
  #qqcomp(fp)
  #quantile(fp, probs = c(0.05, 0.1, 0.2))
  mygof= fitdistrplus::gofstat(fp, discrete=FALSE) # goodness-of-fit statistics
  mykstest = stats::ks.test(imp.vals$t.series, y="pmygpd", location=gpdLocation, scale=gpdScale, shape=0.00001, dummy=1)
  
  par(cex=1.2, bg="white")
  plot(fp, lwd=2, col="cadetblue3")
  
  mtext(side = 1, text = paste0("Page ", numberofplots, " - Declustered - Thunderstorm - POT-GPD Equivalent. Location: ",
                                round(fp$fix.arg$location, digits = 2), ". Scale: ", round(fp$fix.arg$scale, digits = 2),
                                ". Shape: ", round(fp$fix.arg$shape, digits = 2)), outer = TRUE)
  assign(paste0("myprint", numberofplots), recordPlot())
  saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0(outputpath, "myprint", numberofplots, ".rds"))								
  numberofplots = numberofplots + 1
  
  #RMSE
  mx <- seq(min(imp.vals$t.series), max(imp.vals$t.series), length = 10000)
  empCDF <- stats::ecdf(imp.vals$t.series)
  #D of Kolmogorov Smirnov
  #max(abs(empCDF(mx) - pmygpd(mx, scale = fp$fix.arg$scale, location =fp$fix.arg$location, shape = 0.00001, dummy=1)))
  #rmse
  rmse = sqrt(mean((empCDF(mx) - pmygpd(mx, scale = fp$fix.arg$scale, location =fp$fix.arg$location, shape = 0.00001, dummy=1))^2))
  
  myfitgpd = list("Fit_fitdistrplus::fitdist" = capture.output(fp),
                  "GoodnessOfFit_fitdistrplus::gofstat"=capture.output(mygof),
                  "KolmogorovSmirnovTest_stats::ks.test"=capture.output(mykstest),
                  "RMSE"=capture.output(rmse))
  write.xlsx(capture.output(myfitgpd), file=fnfitted, sheetName="t_POT-GPD-Equivalent", append=TRUE, row.names=TRUE)
  
  ##New Idea - Alexys
  #Probability density with the join density according to NIST.SP.500-301.pdf equation 2, page 15.
  #Same equation in Smith 2004. Extreme Values in Finance, Telecommunications, and the Environment,
  #Equation 1.18
  
  #Ip <- function(x) {1/z4 *(1+0.00001*((x-z3)/z4))^((-1/0.00001)-1)}
  
  #joindensityPoissonProcess <- function (serie, serie.dt){  #y1 and y2 are velocities (min and max in dataset),
  #                                                        #t1 and t2 are times (min and max in dataset)
  #  y1 = min(serie)
  #  y2 = max(serie)
  #  t1 = min(serie.dt)
  #  t2 = max(serie.dt)
  #  intensityvalues = lapply (serie, Ip)
  #  product = prod(intensityvalues)
  #  integrand <- Ip
  #  integrated <- integrate(integrand, lower=y1, upper=y2)$value
  #  joindensity <- product * exp(integrated*-1)
  #    Ip.t(lower)/integrated  # Equation after equation 2, page 15
  #}
  
  #Cumulative distribution function according to NIST.SP.500-301.pdf. Equation 5, page 17
  distributionPoissonProcess <- function(velocity, threshold){  #lower is the value of velocity to integrate with it up to infinite
    integrand <- Ip.t
    integrated_tresh_to_vel <- integrate(integrand, lower=threshold, upper=velocity)$value
    integrated_tresh_to_inf <- integrate(integrand, lower=threshold, upper=Inf)$value
    #distribution <- 1- (integrated_tresh_to_vel/integrated_tresh_to_inf)#NIST.SP.500-301.pdf. Equation 5, page 17
    #Be ware I am not using 1 - () as is stated in equation 5
    distribution <- (integrated_tresh_to_vel/integrated_tresh_to_inf)#NIST.SP.500-301.pdf. Equation 5, page 17
  }
  
  #Calculate cumulative distribution function values for velocities 1 to 600 using Poisson Process
  #DfP = sapply(yvels, distributionPoissonProcess, upper=Inf, threshold=z2)
  DfP = sapply(yv, distributionPoissonProcess, threshold=z2)
  library(RColorBrewer)
  cols <- brewer.pal(9,"Set1")
  #Plot: Page 3 - Cumulative distribution function with intensity function of Poisson Process
  #Plot Dfp- Cumulative distribution function  with intensity function of Poisson Process
  plot(x= yv, y=DfP, xlab="Declustered - Thunderstorm - Velocities [Km/h]", ylab="Distribution Function - Poisson Process",
       main=paste("Cumulative Distribution Function from Intensity Function of Poisson Process\n", "Location=", round(z3,2), " Scale=", round(z4,2), "\n", "Station:", number),
       type="l", col=cols[3], lwd=4)
  mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  #h = hist(imp.vals$t.series, probability = TRUE, plot=FALSE)
  #h$counts <- cumsum(h$counts)/sum(h$counts)
  #plot(h, add=TRUE)
  plot(stats::ecdf(imp.vals$t.series), col="cadetblue3", add=TRUE)
  lines(x= yv, y=DfP, col=cols[3], lwd=4)
  curve(pgpd_pp(x, location=z3, shape=0.00000001, scale=z4, threshold=z2), col=cols[5], add=TRUE, lwd=2, lty= 2)
  legend("bottomright", c("Data-Empirical", "Fitted-Theoretical POT-Poisson Process", "Fitted-Theoretical POT-GPD Equivalence"),
         bty = "o",
         col = c("cadetblue3", cols[3], cols[5]),
         lty = c(1, 1, 2), lwd = c(1, 4, 2),
         pch = c(21, NA, NA),
         pt.bg = c("cadetblue3", NA, NA),
         pt.cex = 1, box.lwd = 0, box.col = "white", bg = "white")
  box(lty = 1, col = 'black', lwd=0.5)
  assign(paste0("myprint", numberofplots), recordPlot())
  saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0(outputpath, "myprint", numberofplots, ".rds"))  
  numberofplots = numberofplots + 1
  
  #Plot: Page 4
  #Histogram and fitted density function from density function of NIST.SP.500-301.pdf. Equation after equation 2, page 15
  #Whitout Substracting Threshold
  overthreshold = imp.vals$t.series >= z2 #Be sure only data above threshold are included
  #hist(imp.vals$t.series[overthreshold], prob=TRUE, main=paste("Fitted density function from Intensity Function of Poisson Process\nWhitout Substracting Threshold\nStation:", number, sep=" "))
  #lines(x=yvels, y=dfP, col="red")
  #mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  #numberofplots = numberofplots + 1
  
  #Plot: Page 5
  #Histogram and fitted density function from density function of NIST.SP.500-301.pdf. Equation after equation 2, page 15
  #Subtracting Threshold
  #overthreshold = imp.vals$t.series >= z2 #Be sure only data above threshold are included
  #hist(imp.vals$t.series[overthreshold]-nt.thresh, prob=TRUE, main=paste("Fitted density function from Intensity Function of Poisson Process\nSubstracting Threshold\nStation:", number, sep=" "))
  #lines(x=yvels, y=dfP, col="red")
  #mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  #numberofplots = numberofplots + 1
  
  #Calculate velocities for typical return periods using Poisson Process
  veocitiesfortypicalreturnperiodsP <- approx(x=paP, y=yvels, xout = typicalExcedenceProbabilities)$y  #Interpolate tipical excedence probabilities to get velocities using Poisson
  
  zzz[zz,10:20] = veocitiesfortypicalreturnperiodsP #Velocities for typical return periods Poisson
  
  #tipicalReturnPeriods = c("10","20","50","100","250","500","700","1000","1700","3000","7000")
  #
  #Plot: Page 6
  maxy = max(veocitiesfortypicalreturnperiodsP) + 0.15*max(veocitiesfortypicalreturnperiodsP)
  #Hazards Curves using Intensity Function of Poisson Process
  library(RColorBrewer)
  cols <- brewer.pal(9,"Set1")
  
  plot(x= 1/paP, y=yvels, xlab="Return Periods (Years) - Poisson Process Intensity Function",
       ylab="Velocities Km/h", main=paste("Declustered - Thunderstorms - Hazard Curve\n", "Location=", round(z3,2), " Scale=", round(z4,2), "\n", "Station:", number),
       xlim=c(0,10000), ylim = c(0, maxy), type="l", col=cols[3], lwd=3)
  #curve (returnlevelPP(x, location=z3, shape=0.00001, scale=z4), col="red", add=TRUE, lwd=6)
  myx= 1:10000
  myrlpp = sapply(myx, returnlevelPP, location=z3, shape=0.00001, scale=z4)
  lines(myx, myrlpp, col=cols[4], lwd=2, lty=5)
  
  #curve (returnlevelGPD(x, location=z3, shape=0.00001, scale=z4, threshold=z2), col="black", add=TRUE, lwd=2)
  myrlgpd2 = sapply(myx, returnlevelGPD, location=z3, shape=0.00001, scale=z4, threshold=z2)
  lines(myx, myrlgpd2, col=cols[5], lwd=1, lty=2)
  points(x= tipicalReturnPeriods, y= veocitiesfortypicalreturnperiodsP, col=cols[9], pch=20)
  text(x = tipicalReturnPeriods, y = veocitiesfortypicalreturnperiodsP,
       labels = paste0("(",tipicalReturnPeriods,",",round(veocitiesfortypicalreturnperiodsP, digits=1),")"),
       cex=0.8, pos = 4)
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
  saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0(outputpath, "myprint", numberofplots, ".rds"))
  numberofplots = numberofplots + 1
  
  
  #rlgpd = sapply(1:7000, returnlevelGPD, location=z3, shape=0.00001, scale=z4, threshold=z2)
  #plot(x= 1:7000, y=rlgpd, col="black", type="l")
  #plot(x= myx, y=myrlgpd2, col="black", type="l")
  
  #___
  #___b) Poisson Process - Gumbel like tile length. Shape = 0
  #___
  
  
  #Intensity Function Poisson Process - Gumbel like tile length. Shape = 0. NIST.SP.500-301.pdf, page 32, equation 10
  Ig.t <- function(x) {1/z4 * exp(-(x-z3)/z4)}
  at = imp.vals$n.thunders.per.year*(1/24)
  #Representing equation 4, page 16 (NIST.SP.500-301.pdf), but using Ip
  returnlevel_equation_gumbel <- function(lower, upper, at){  #lower is the value of velocity to integrate with it up to infinite
    integrated.t <- integrate(Ig.t, lower=lower, upper=upper)$value
    value = at * integrated.t
  }
  
  #Calculate excedence probability for velocities 1 to 600 using Poisson Process - Gumbel Intensity Function
  paG = sapply(yvels, returnlevel_equation_gumbel, upper=Inf, at=at)
  
  #Calculate velocities for typical return periods using Poisson Process - Gumbel Intensity Function
  veocitiesfortypicalreturnperiodsG <- approx(x=paG, y=yvels, xout = typicalExcedenceProbabilities)$y  #Interpolate tipical excedence probabilities to get velocities using Gumbel
  zzz[zz,21:31] = veocitiesfortypicalreturnperiodsG #Velocities for typical return periods Gumbel Tile Length Intensity Function
  
  #Plot: Page 7
  #Hazards Curves using Intensity Function of Gumbel like tail length
  maxy = max(veocitiesfortypicalreturnperiodsG) + 0.15*max(veocitiesfortypicalreturnperiodsG)
  library(RColorBrewer)
  cols <- brewer.pal(9,"Set1")
  plot(x= 1/paG, y=yvels, xlab="Return Periods (Years) - Gumbel like tail Intensity Function of Poisson Process",
       ylab="Velocities Km/h", main=paste("Declustered - Thunderstorms - Hazard Curve\n", "Location=", round(z3,2), " Scale=", round(z4,2), "\n", "Station:", number),
       xlim=c(0,10000), ylim=c(0, maxy), type="l", col=cols[3], lwd=3)
  myx= 1:10000
  myrlpp = sapply(myx, returnlevelPP, location=z3, shape=0.00001, scale=z4)
  lines(myx, myrlpp, col=cols[4], lwd=2, lty=5)
  #curve (returnlevelPP(x, location=z3, shape=0.00001, scale=z4), col="red", add=TRUE, lwd=6)
  myrlgpd2 = sapply(myx, returnlevelGPD, location=z3, shape=0.00001, scale=z4, threshold=z2)
  lines(myx, myrlgpd2, col=cols[5], lty=2, lwd=1)
  #curve (returnlevelGPD(x, location=z3, shape=0.00001, scale=z4, threshold=z2), col="black", add=TRUE, lwd=2)
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
  saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0(outputpath, "myprint", numberofplots, ".rds"))  
  numberofplots = numberofplots + 1
  
  
  #___
  #___c) Using Gumbel Quantile replacing fitted parameters (location=z3, scale=z4)
  #___
  
  
  
  #Probability density function using Gumbel distribution
  #Plot: Page 8
  
  library(RcmdrMisc)
  .x <- seq(-300, 300, length.out=1000)
  dfG = RcmdrMisc::dgumbel(.x, location=z3, scale=z4)
  plot(.x, dfG, col="red",
       xlab="Thunderstorms - Velocities Km/h", ylab="Density Function - Gumbel Distribution",
       main=paste("Gumbel Density Function, but using parameters of Poisson Process\n", "Location=", round(z3,2), " Scale=", round(z4,2), "\n", "Station:", number),
       type="l", lwd=1, ylim=c(0,0.08))
  #plotDistr(.x, dgumbel(.x, location=0, scale=1), cdf=FALSE, xlab="x",
  #          ylab="Density", main=paste("Gumbel Distribution:  Location=0, Scale=1"))
  mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  
  pdfG <- function(x) {1/z4 *exp(-(x-z3)/z4)*exp(-exp(-(x-z3)/z4))}
  #curve(pdfG, add=T)
  hist(imp.vals$t.series[overthreshold], prob=TRUE, add=T, col="cadetblue3")
  lines(.x, dfG, col="red")
  legend("topright", c("Data-Empirical", "Fitted-Theoretical"),
         bty = "n",
         col = c("cadetblue3", "red"),
         lty = c(0, 1), lwd = c(0, 1),
         pch = c(22, NA),
         pt.bg = c("cadetblue3", NA),
         pt.cex = 2)
  assign(paste0("myprint", numberofplots), recordPlot())
  saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0(outputpath, "myprint", numberofplots, ".rds"))		 
  numberofplots = numberofplots + 1
  
  #pdfG <- function(x) {1/z4 *exp(-(x-0)/z4)*exp(-exp(-(x-0)/z4))}
  #curve(pdfG, add=T)
  
  
  #Cumulative distribution function with Gumbel distribution
  #Plot: Page 9
  .x <- seq(-300, 300, length.out=1000)
  DfG = RcmdrMisc::pgumbel(.x, loc=z3, scale=z4)
  plot(.x, DfG,
       xlab="Thunderstorms - Velocities Km/h", ylab="Cumulative Distribution Function - Gumbel Distribution",
       main=paste("Gumbel Cumulative Distribution, but using parameters of Poisson Process\n", "Location=", round(z3,2), " Scale=", round(z4,2), "\n", "Station:", number),
       type="l", col="red", lwd=2)
  #lines(yvels, DfG, col="red", lty=4)
  mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  plot(stats::ecdf(imp.vals$t.series), col="cadetblue3", add=TRUE)
  lines(.x, DfG, col="red")
  legend("bottomright", c("Data-Empirical", "Fitted-Theoretical"),
         bty = "o",
         col = c("cadetblue3", "red"),
         lty = c(1, 1), lwd = c(1, 2),
         pch = c(21, NA),
         pt.bg = c("cadetblue3", NA),
         pt.cex = 1, box.lwd = 0, box.col = "white", bg = "white")
  box(lty = 1, col = 'black', lwd=0.5)
  assign(paste0("myprint", numberofplots), recordPlot())
  saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0(outputpath, "myprint", numberofplots, ".rds")) 
  numberofplots = numberofplots + 1
  
  #Plot: Page 10
  #Histogram and fitted density function from Gumbel distribution
  overthreshold = imp.vals$t.series >= z2 #Be sure only data above threshold are included
  #imp.vals$t.series[imp.vals$t.series > t.thresh]
  hist(imp.vals$t.series[overthreshold], prob=TRUE, xlab="Declustered Thunderstorm", col="cadetblue3",
       main=paste("Fitted Gumbel density function using parameters of Poisson Process\n", "Location=", round(z3,2), " Scale=", round(z4,2), "\n", "Station:", number))
  curve(RcmdrMisc::dgumbel(x, location=z3, scale=z4), add=TRUE, col="red", lwd=2)
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
  saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0(outputpath, "myprint", numberofplots, ".rds"))		 
  numberofplots = numberofplots + 1
  
  #Following section "Using the Fitted Model" of Extreme Wind Speeds: Overview (https://www.itl.nist.gov/div898/winds/overview.htm)
  #Calculation of Velocities for Return Periods using Percent Point Function of Gumbel (Quantiles: qgumbel)
  
  gumbelVelocitiesQuantileFunction <- function(mri){  #mri: return intervals from 1 to 3000
    RcmdrMisc::qgumbel((1-(1/mri)), location=z3, scale=z4)
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
    RcmdrMisc::qgumbel((1-(1/(lambda*mri))), location=z3, scale=z4)
  }
  
  gumbelVelocitiesQuantileFunctionNpy <- function(mri, npy){  #mri: return intervals from 1 to 3000
    RcmdrMisc::qgumbel((1-(1/(npy*mri))), location=z3, scale=z4)
  }
  
  
  gumbelVelocitiesQuantileFunctionBerry <- function(mri, npy, truncate){  #mri: return intervals from 1 to 3000
    probs = 1 -(1/(mri*npy))
    probs2 = (probs-truncate)/(1-truncate)
    RcmdrMisc::qgumbel(probs2, location=z3, scale=z4)
  }
  
  #Calculation of Velocities for Typical Return Periods using Quantile Function of Gumbel Distribution
  #veocitiesfortypicalreturnperiodsQG = sapply(tipicalReturnPeriods, function(x){qgumbel((1-(1/x)), location=z3, scale=z4)})
  veocitiesfortypicalreturnperiodsQG = sapply(tipicalReturnPeriods, gumbelVelocitiesQuantileFunction)
  
  #Next option whitout using Lambda
  #npy=imp.vals$n.thunders.per.year
  #veocitiesfortypicalreturnperiodsQG = sapply(tipicalReturnPeriods, gumbelVelocitiesQuantileFunctionLambda,
  #                                            npy=npy, numberofsamples= length(imp.vals$t.series),
  #                                            numberofsamplesoverthreshold= length(imp.vals$t.series[overthreshold]) )
  
  #This option using npy
  #npy=imp.vals$n.thunders.per.year
  #veocitiesfortypicalreturnperiodsQG = sapply(tipicalReturnPeriods, gumbelVelocitiesQuantileFunctionNpy, npy=npy)
  
  
  #This option using berry
  #npy=imp.vals$n.thunders.per.year
  #truncate = 1 - (sum(overthreshold)/length(imp.vals$t.series))
  #veocitiesfortypicalreturnperiodsQG = sapply(tipicalReturnPeriods, gumbelVelocitiesQuantileFunctionBerry, npy=npy, truncate=truncate)
  
  
  
  zzz[zz,32:42] = veocitiesfortypicalreturnperiodsQG #Velocities for typical return periods using Gumbel quantile function
  
  
  #Calculation of velocities for return periods from 1 to 3000
  mri = 1:10000
  allVelGumbelQuantileFunction = sapply(mri, gumbelVelocitiesQuantileFunction)
  library(RColorBrewer)
  cols <- brewer.pal(9,"Set1")
  #Hazard curve using quantile function of Gumbel distribution
  #Plot: Page 11
  plot(x= mri, y=allVelGumbelQuantileFunction, xlab="Return Periods (Years) - Gumbel Quantile Function using parameters of Poisson Process",
       ylab="Velocities Km/h", main=paste("Declustered - Thunderstorms - Hazard Curve\n", "Location=", round(z3,2), " Scale=", round(z4,2), "\n", "Station:", number),
       xlim=c(0,10000), type="l", col=cols[3], lwd=3)
  myx= 1:10000
  myrlpp = sapply(myx, returnlevelPP, location=z3, shape=0.00001, scale=z4)
  lines(myx, myrlpp, col=cols[4], lwd=2, lty=5)
  #curve (returnlevelPP(x, location=z3, shape=0.00001, scale=z4), col="red", add=TRUE, lwd=6)
  myrlgpd2 = sapply(myx, returnlevelGPD, location=z3, shape=0.00001, scale=z4, threshold=z2)
  lines(myx, myrlgpd2, col=cols[5], lwd=1, lty=2)
  #curve (returnlevelGPD(x, location=z3, shape=0.00001, scale=z4, threshold=z2), col="black", add=TRUE, lwd=2)
  points(x= tipicalReturnPeriods, y= veocitiesfortypicalreturnperiodsQG, col=cols[9], pch=20)
  text(x = tipicalReturnPeriods, y = veocitiesfortypicalreturnperiodsQG, labels = paste0("(",tipicalReturnPeriods,",",round(veocitiesfortypicalreturnperiodsQG, digits=1),")"), cex=0.8, pos = 4)
  mtext(side = 1, text = paste0("Page ", numberofplots), outer = TRUE)
  legend("bottomright", c("Gumbel Quantile","POT-PP (formula)", "POT-GPD equivalence"),
         bty = "n",
         col = cols[3:5],
         lty = c(1, 5, 2), lwd = c(3, 2, 1),
         pch = c(NA, NA, NA),
         pt.bg = c(NA, NA, NA),
         pt.cex = c(0,0,0))
  box(lty = 1, col = 'black', lwd=0.5)
  assign(paste0("myprint", numberofplots), recordPlot())
  saveRDS(eval(parse(text=paste0("myprint", numberofplots))), paste0(outputpath, "myprint", numberofplots, ".rds"))  
  numberofplots = numberofplots + 1
}
