### RD design functions
###peter crosta pmcrosta at gmail dot com
## This file consists of only functions that are used in the 
## RD analysis. It will be sourced

###########MOST FUNCTIONS ARE DEPENDENT ON maxbin AND binwidth (h) ###########################

## Create class binlist for numeric breaks
setClass("binlist", representation(breaks="numeric", right="numeric", left="numeric"))

getbins <- function(binwidth, maxbin) {
  rightbins <- seq(binwidth, maxbin, binwidth)
  leftbins <- -rev(seq(0, maxbin, binwidth))

  breakbins <- c(sort(leftbins), rightbins)
  
  z <- new("binlist")
  z@breaks <- breakbins[abs(breakbins)<=maxbin]
  z@right <- rightbins[rightbins<=maxbin]
  z@left <- leftbins[leftbins>=-maxbin]

  return(z)
}

# -244 to 121 might matter, days before and after Sep 1
mainhist <- function(binwidth, maxbin, yesplot=TRUE, emcp=emc) {
  breaks <- getbins(binwidth, maxbin)@breaks
  emcdum <- ifelse(emcp >= min(breaks) & emcp <= max(breaks), 1, 0)
  maxbin.hist <- hist(emcp[emcdum==1], breaks=breaks, plot=yesplot, main=paste("Histogram of birth date (centered): Binwidth=",binwidth, ", ", length(breaks)-1, " bins", sep=""), xlim=c(min(breaks), max(breaks)), xlab="Centered Birth Date")
  if (yesplot) abline(v=0, lty=2, lwd=2, col="darkgrey")
  binmids <- maxbin.hist$mids  
}


####this function is built from McCrary(2008); tests discontinuity of density
mccrary <- function(binwidth, maxbin, nn, emcp=emc) {
    f.triangle <- function(x) { ifelse(abs(x)<1, 1-abs(x), 0) }

    breaks <- getbins(binwidth, maxbin)@breaks
    emcdum <- ifelse(emcp >= min(breaks) & emcp <= max(breaks), 1, 0)
    maxbin.hist <- hist(emcp[emcdum==1], breaks=breaks, plot=FALSE)
    binmids <- maxbin.hist$mids
    
    totaln <- sum(maxbin.hist$counts)
    Xj <- binmids
    Yj <- maxbin.hist$counts / (totaln*binwidth)

    ##set bandwidth (nearest neighbors); #nn <-20; this is argument in func call

    ##first step histogram
    plot(Xj, Yj, xlab="Bin midpoints", ylab="Normalized cell size", main=paste("First step histogram for McCrary(2008)\nbinwidth=",binwidth, " bandwidth=", nn, sep=""))
    
    ##local fit right and plot
    right.locfit <- locfit.raw(lp(Xj[Xj > 0], h=nn, deg=1), Yj[Xj > 0], kern='tria', family="gaussian")
    plot(right.locfit, add=TRUE, col="red")
    
    ###weighted reg around boundary point with triang kernel to get intercept term 
    wt <- f.triangle((Xj[Xj > 0]-rep(0, length( Xj[Xj > 0])))/nn)
    wt <- (wt/(sum(wt)))*length(wt)
    templm <- lm(Yj[Xj > 0]~Xj[Xj > 0], weight=wt)
    fr <- coef(templm)[1]
    
    #local fit left and plot
    left.locfit <- locfit.raw(lp(Xj[Xj <= 0], h=nn, deg=1), Yj[Xj <= 0], kern='tria', family="gaussian")
    plot(left.locfit, add=TRUE, col="red")
    
    ###weighted reg around boundary point with trian kernel to get intercept term
    wt <- f.triangle((Xj[Xj <= 0]-rep(0,length(Xj[Xj <= 0])))/nn)
    wt <- (wt/(sum(wt)))*length(wt) #like aweights in stata
    templm <- lm(Yj[Xj <= 0]~Xj[Xj <= 0], weight=wt)
    fl <- coef(templm)[1]
    
    abline(v=0, lty=2, lwd=2, col="darkgrey")
    
    thetahat <- log(fr) - log(fl)
    sethetahat <- sqrt( (1/(totaln*nn)) * (24/5) * ((1/fr) + (1/fl)) )

    legend("topright", bg="white", cex=.9, pt.cex=.9, c("Cutoff", "Bin-mid regression with triangular kernel"), lty=c(2, 1), lwd=c(2, 1), col=c("darkgrey", "red"))
    text(-maxbin+30, max(Yj)*.92, paste("theta-hat =", round(thetahat, 4), "\nstd. err. = ", round(sethetahat,4)))
    return(cat(binwidth, "\t", maxbin, "\t", nn, "\t", thetahat, "\t", sethetahat, "\n"))
}


#####function to generate prop tests to look at subgroup covariate differences around cutoff (like gifted and atrisk)
prettyprop <- function(binwidth, maxbin, covariate, emcp=emc, lab="") {
    breaks <- getbins(binwidth, maxbin)@breaks
    emcdum <- ifelse(emcp >= min(breaks) & emcp <= max(breaks), 1, 0)
 
    covariate <- covariate[emcdum==1]
    cutbybreak <- cut(emcp[emcdum==1], breaks=breaks)

    afters<-which(cutbybreak %in% tail(levels(cutbybreak), maxbin))
    befores<-which(cutbybreak %in% head(levels(cutbybreak), maxbin))
    
    x1 <- sum(covariate[befores], na.rm=TRUE)
    x2 <- sum(covariate[afters], na.rm=TRUE)
    n1 <- length(na.omit(covariate[befores]))
    n2 <- length(na.omit(covariate[afters]))
  
    pt <- prop.test(c(x1,x2), c(n1, n2))
    return(paste(lab, round(pt$estimate[1], 3), round(pt$estimate[2], 3), round(diff(pt$estimate), 3), round(pt$p.value, 3), n1, n2, sep=', '))
}

######## similar to prettyprop, but for variables that take a t-test
ttest <- function(binwidth, maxbin, covariate, emcp=emc, lab="") {  
    breaks <- getbins(binwidth, maxbin)@breaks
    emcdum <- ifelse(emcp >= min(breaks) & emcp <= max(breaks), 1, 0)
 
    covariate <- covariate[emcdum==1]
    cutbybreak <- cut(emcp[emcdum==1], breaks=breaks)
    
    afters<-which(cutbybreak %in% tail(levels(cutbybreak), maxbin))
    befores<-which(cutbybreak %in% head(levels(cutbybreak), maxbin))
    
    x1 <- covariate[befores]
    x2 <- covariate[afters]
    tt <- t.test(x1, x2)

    return(paste(lab, round(tt$estimate[1], 3), round(tt$estimate[2], 3), round(-diff(tt$estimate), 3), round(tt$p.value, 3), length(na.omit(x1)), length(na.omit(x2)), sep=","))
}


covfunc <- function(binwidth, maxbin, covariate, ylabcov, use.loess=FALSE, emcp=emc, YLIM=TRUE, nn=10, g) {
  
    breaks <- getbins(binwidth, maxbin)@breaks
    emcdum <- ifelse(emcp >= min(breaks) & emcp <= max(breaks), 1, 0)
    
    cutbybreak <- cut(emcp[emcdum==1], breaks=breaks)
    binmids <- hist(emcp[emcdum==1], breaks=breaks, plot=FALSE)$mids
    
    binned.covariate <- tapply(covariate[emcdum==1], cutbybreak, function(x) mean(x, na.rm=TRUE))
    forcing <- ifelse(binmids>0, 1, 0)
    tx <- which(forcing==1)
    forcing.sample <- ifelse(emcp>0, 1, 0)
    
    if (YLIM) {   
        plot(binned.covariate~binmids, col=NULL, xlab="Birth date (centered)", ylab=ylabcov, main=paste(ylabcov, " Grade ",g, "\nBinwidth=",binwidth, ", ", length(breaks)-1, " bins", sep=""), xaxt='n', ylim=c(0, 1))
    } else plot(binned.covariate~binmids, col=NULL, xlab="Birth date (centered)", ylab=ylabcov, main=paste(ylabcov, " Grade ",g, "\nBinwidth=",binwidth, ", ", length(breaks)-1, " bins", sep=""), axes=FALSE)
    
    points(binmids[tx], binned.covariate[tx], col="blue4", pch=4)
    points(binmids[-tx], binned.covariate[-tx], col="darkorange", pch=16)
    axis(1, at=seq(-maxbin, maxbin, as.integer(maxbin/10)))
    
    if (!YLIM) axis(2, at=round(seq(min(binned.covariate, na.rm=TRUE), max(binned.covariate, na.rm=TRUE), diff(range(binned.covariate, na.rm=TRUE))/10),2))
    
    abline(v=0, lty=2, lwd=2, col="darkgrey")
    
    if (ylabcov %in% covstlab) {
      pval <- ttest(1, binwidth, covariate, emcp)
    } else {
      pval <- prettyprop(1, binwidth, covariate, emcp)
    }
    pvalsplit <- as.numeric(unlist(strsplit(pval, split=",")))
    cat(pvalsplit, "\n")
    
    if (use.loess) {      
      f.triangle <- function(x) { ifelse(abs(x)<1, 1-abs(x), 0) }
      Yj <- binned.covariate
      Xj <- binmids[which(!is.na(Yj))]
      Yj <- binned.covariate[which(!is.na(Yj))]      
      
      #changed so this is now argument
      #nn<-binwidth*3
        
      ##local fit right and plot
      right.locfit <- locfit.raw(lp(Xj[Xj > 0], h=nn, deg=1), Yj[Xj > 0], kern="tria", family="gaussian")
      plot(right.locfit, add=TRUE, col="green")    
      #local fit left and plot
      left.locfit <- locfit.raw(lp(Xj[Xj <= 0], h=nn, deg=1), Yj[Xj <= 0], kern='tria', family="gaussian")
      plot(left.locfit, add=TRUE, col="green")   
      
      text(-maxbin/2, quantile(binned.covariate, na.rm=TRUE)[4], cex=0.8, paste("1 binwidth difference at cutoff =", round(pvalsplit[4],3), "\np-val=", round(pvalsplit[5],3)))
      #Add legend.
      legend("topright", bg="white", cex=.9, pt.cex=.9, c("Older", "Younger", "Cutoff", "Piecewise loess"), lty=c(NA, NA, 2, 1), lwd=c(NA, NA, 2, 1), col=c("blue4", "darkorange", "darkgrey", "green"), pch=c(4, 16, NA, NA, NA))    
    }  else {
        ##run simple regression model and add curves to plot
        model <- lm(covariate~forcing.sample + emcp + I(emcp^2), subset=emcdum==1)
        coefs <- coefficients(model)
        curve(coefs[1] + coefs[3]*x + coefs[4]*x^2, -maxbin, 0, add=T, col="red")
        curve(coefs[1] + coefs[2] + coefs[3]*x + coefs[4]*x, 0, maxbin, add=T, col="red")
        legend("topright", bg="white", cex=.9, pt.cex=.9, c("Older", "Younger", "Cutoff", "Piecewise OLS"), lty=c(NA, NA, 2, 1), lwd=c(NA, NA, 2, 1), col=c("blue4", "darkorange", "darkgrey", "red"), pch=c(4, 16, NA, NA))
        text(-maxbin/2, quantile(binned.covariate, na.rm=TRUE)[4], cex=0.8, paste("OLS difference at cutoff = ", round(coefs[2], 3), "(", round(sqrt(vcov(model)[2,2]), 3), ")",sep=""))    
    }
    
}


outvforce <- function(binwidth, maxbin, covariate, ylabcov, use.loess=FALSE, emcp=emc, YLIM=TRUE, nn=10, g, change=FALSE) {  

  ###Plot of outcome vs. forcing variable
  breaks <- getbins(binwidth, maxbin)@breaks
  emcdum <- ifelse(emcp >= min(breaks) & emcp <= max(breaks), 1, 0)

  cutbybreak <- cut(emcp[emcdum==1], breaks=breaks)
  binmids <- hist(emcp[emcdum==1], breaks=breaks, plot=FALSE)$mids

  ##group means with binwidth of ; midpoints for horiz axis supplied by hist
  binned.covariate <- tapply(covariate[emcdum==1], cutbybreak, function(x) mean(x, na.rm=TRUE))
  forcing <- ifelse(binmids>0, 1, 0)
  tx <- which(forcing==1)
  forcing.sample <- ifelse(emcp>0, 1, 0)

  if (change) {
    maintit=paste(ylabcov, " Grade ",g-1, " to ", g, "\nBinwidth=",binwidth, ", ", length(breaks)-1, " bins", sep="")
  } else maintit=paste(ylabcov, " Grade ",g, "\nBinwidth=",binwidth, ", ", length(breaks)-1, " bins", sep="")
  
  if (YLIM) {
        plot(binned.covariate~binmids, col=NULL, xlab="Birth date (centered)", ylab=ylabcov, main=maintit, xaxt='n', ylim=c(0, 1))
  } else plot(binned.covariate~binmids, col=NULL, xlab="Birth date (centered)", ylab=ylabcov, main=maintit, axes=FALSE)

  points(binmids[tx], binned.covariate[tx], col="blue4", pch=4)
  points(binmids[-tx], binned.covariate[-tx], col="darkorange", pch=16)
  axis(1, at=seq(-maxbin, maxbin, as.integer(maxbin/10)))

  if (!YLIM) axis(2, at=round(seq(min(binned.covariate, na.rm=TRUE), max(binned.covariate, na.rm=TRUE), diff(range(binned.covariate, na.rm=TRUE))/10),2))

  abline(v=0, lty=2, lwd=2, col="darkgrey")


 if (ylabcov %in% covstlab) {
      pval <- ttest(1, binwidth, covariate, emcp)
    } else {
      pval <- prettyprop(1, binwidth, covariate, emcp)
    }
    pvalsplit <- as.numeric(unlist(strsplit(pval, split=",")))
    cat(pvalsplit, "\n")

    ### LAST MOD IS HERE
    if (use.loess) {      
    f.triangle <- function(x) { ifelse(abs(x)<1, 1-abs(x), 0) }
    Yj <- binned.covariate
    Xj <- binmids[which(!is.na(Yj))]
    Yj <- binned.covariate[which(!is.na(Yj))] 

    #changed so this is now argument
    #nn<-binwidth*3    
         
    ##local fit right and plot
    right.locfit <- locfit.raw(lp(Xj[Xj > 0], h=nn, deg=1), Yj[Xj > 0], kern='tria', family="gaussian")
    plot(right.locfit, add=TRUE, col="green")
       
    #local fit left and plot
    left.locfit <- locfit.raw(lp(Xj[Xj <= 0], h=nn, deg=1), Yj[Xj <= 0], kern='tria', family="gaussian")
    plot(left.locfit, add=TRUE, col="green")   
    
    text(-maxbin/2, quantile(binned.covariate, na.rm=TRUE)[4], paste("1 binwidth difference at cutoff =", round(pvalsplit[4],3), "\np-val=", round(pvalsplit[5],3))) 
    #Add legend.
      legend("topright", bg="white", cex=.9, pt.cex=.9, c("Older", "Younger", "Cutoff", "Piecewise loess"), lty=c(NA, NA, 2, 1), lwd=c(NA, NA, 2, 1), col=c("blue4", "darkorange", "darkgrey", "green"), pch=c(4, 16, NA, NA, NA))    
    }  else {
    ##run simple regression model and add curves to plot
    model <- lm(binned.covariate~forcing + binmids + I(binmids^2), subset=emcdum==1)
    #model <- lm(covariate~forcing.sample + emcp + I(emcp^2), subset=emcdum==1)
    coefs <- coefficients(model)
    curve(coefs[1] + coefs[3]*x + coefs[4]*x^2, -maxbin, 0, add=T, col="red")
    curve(coefs[1] + coefs[2] + coefs[3]*x + coefs[4]*x, 0, maxbin, add=T, col="red")
    legend("topright", bg="white", cex=.9, pt.cex=.9, c("Older", "Younger", "Cutoff", "Piecewise OLS"), lty=c(NA, NA, 2, 1), lwd=c(NA, NA, 2, 1), col=c("blue4", "darkorange", "darkgrey", "red"), pch=c(4, 16, NA, NA))
        text(-maxbin/2, quantile(binned.covariate, na.rm=TRUE)[4], cex=0.8, paste("OLS difference\nat cutoff = ", round(coefs[2], 3), "(", round(sqrt(vcov(model)[2,2]), 3), ")",sep=""))    
    }
}

##********************************************************************COVARIATES******************
#x11(height=10, width=10)

clx <- function(fm, dfcw, cluster){
         # R-codes (www.r-project.org) for computing
         # clustered-standard errors. Mahmood Arai, Jan 26, 2008.
   
	 # The arguments of the function are:
         # fitted model, cluster1 and cluster2
         # You need to install libraries `sandwich' and `lmtest'
	 
         # reweighting the var-cov matrix for the within model
         #library(sandwich);library(lmtest)
         M <- length(unique(cluster))   
         N <- length(cluster)           
         K <- fm$rank                        
         dfc <- (M/(M-1))*((N-1)/(N-K))  
         uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
         vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)*dfcw
         coeftest(fm, vcovCL) 
}

mclx <- function(fm, dfcw, cluster1, cluster2){
         # R-codes (www.r-project.org) for computing multi-way 
         # clustered-standard errors. Mahmood Arai, Jan 26, 2008. 
         # See: Thompson (2006), Cameron, Gelbach and Miller (2006)
         # and Petersen (2006).
	 # reweighting the var-cov matrix for the within model
         
         # The arguments of the function are:
         # fitted model, cluster1 and cluster2
         # You need to install libraries `sandwich' and `lmtest'
         
         library(sandwich);library(lmtest)
         cluster12 = paste(cluster1,cluster2, sep="")
         M1  <- length(unique(cluster1))
         M2  <- length(unique(cluster2))   
         M12 <- length(unique(cluster12))
         N   <- length(cluster1)          
         K   <- fm$rank             
         dfc1  <- (M1/(M1-1))*((N-1)/(N-K))  
         dfc2  <- (M2/(M2-1))*((N-1)/(N-K))  
         dfc12 <- (M12/(M12-1))*((N-1)/(N-K))  
         u1j   <- apply(estfun(fm), 2, function(x) tapply(x, cluster1,  sum)) 
         u2j   <- apply(estfun(fm), 2, function(x) tapply(x, cluster2,  sum)) 
         u12j  <- apply(estfun(fm), 2, function(x) tapply(x, cluster12, sum)) 
         vc1   <-  dfc1*sandwich(fm, meat=crossprod(u1j)/N )
         vc2   <-  dfc2*sandwich(fm, meat=crossprod(u2j)/N )
         vc12  <- dfc12*sandwich(fm, meat=crossprod(u12j)/N)
         vcovMCL <- (vc1 + vc2 - vc12)*dfcw
         coeftest(fm, vcovMCL)
}

leecard <- function(binwidth, maxbin, outcome, testtype, emcp=emc) {
    ## outcome = const + B1*age + B2*Cutoff + e
    ## outcome = const + B1*emcp + B2*forcing.sample + e

  breaks <- getbins(binwidth, maxbin)@breaks
  emcdum <- ifelse(emcp >= min(breaks) & emcp <= max(breaks), 1, 0)

  cutbybreak <- cut(emcp[emcdum==1], breaks=breaks)
  binmids <- hist(emcp[emcdum==1], breaks=breaks, plot=FALSE)$mids

  emcbinmids <- findInterval(emcp[emcdum==1], binmids, rightmost.closed=TRUE)
  emcmids <-  binmids[emcbinmids]
  
  #binned.outcome <- tapply(outcome[emcdum==1], cutbybreak, function(x) mean(x, na.rm=TRUE))
  #binned.variance <- tapply(outcome[emcdum==1], cutbybreak, function(x) var(x, na.rm=TRUE))
  
  #binned.test <- tapply(outcome[emcdum==1], cutbybreak, function(x) table(!is.na(x)))
  #Nj <- binned.takers <- unlist(lapply(binned.test, function(x) sum(x["TRUE"])))
  #binned.takers <- binned.takers / (sum(binned.takers) / length(binned.takers))
  #nomissbin <- sapply(names(unlist(binned.takers)), function(x) substr(x, start=1, stop=nchar(x)-5))
  
  forcing <- ifelse(binmids>0, 1, 0)
  tx <- which(forcing==1)
  forcing.sample <- ifelse(emcp>0, 1, 0)

  #leftcut <- sort(getbins(binwidth, maxbin)@left)[1]
  #rightcut <- sort(getbins(binwidth, maxbin)@right, decreasing=TRUE)[1]
  #rightobs <- ifelse(emc > 0 & emc <= rightcut, 1, 0)
  #leftobs <- ifelse(emc <= 0 & emc > leftcut, 1, 0)
  #obs <- ifelse(emc >= leftcut & emc <= rightcut, 1, 0)
  test <- outcome[emcdum==1][emcbinmids!=0]
  alpha <- forcing.sample[emcdum==1][emcbinmids!=0]

  ##for GOF test
  Unrestricted <- lm(test~alpha+as.factor(emcmids))
  #coef on forcing.sample is NA
  J <- Unrestricted$rank
  N <- dim(Unrestricted$model)[1]
  ESSur <- sum(Unrestricted$residuals^2)
  
  
  ####global polynomial models over microdata with cluster consistent standard errors and a GOF test (Card and Lee 2008)
  #p0.T <- lm(outcome~forcing.sample, subset=emcdum==1)
  #clx(p0.T, 1, clustersub)  #cluster consistent standard errors
  #sqrt(vcovHC(p0.T, "HC")) #robust Heteroskedasticity-consistent standard errors
  #K <- p0.T$rank
  #ESSr <- sum(p0.T$residuals^2)
  #G = ((ESSr-ESSur)/(J-K)) / (ESSur / (N-J)) # ~F(J-K, N-J)
  #print(pval <- pf(G, J-K, N-J))  GOODNESS OF FIT TEST 
  
  cat("LeeCard GOF Tests", "P-values from F-test", "Polynomial order\n", sep=", ")
  p1.T <- lm(test~alpha+emcmids)
  clustersub <- emcp[attributes(p1.T$model)$row.names]
  c1<-clx(p1.T, 1, clustersub)
  K <- p1.T$rank 
  ESSr <- sum(p1.T$residuals^2)
  G = ((ESSr-ESSur)/(J-K)) / (ESSur / (N-J)) # ~F(J-K, N-J)
  cat(pval <- pf(G, J-K, N-J), K-2, "\n", sep=", ") #
    
  p2.T <- p3.T <- update(p1.T, . ~ . + I(emcmids^2))
  clustersub <- emcp[attributes(p2.T$model)$row.names]
  c2<-clx(p2.T, 1, clustersub)
  K <- p2.T$rank 
  ESSr <- sum(p2.T$residuals^2)
  G = ((ESSr-ESSur)/(J-K)) / (ESSur / (N-J)) # ~F(J-K, N-J)
  cat(pval <- pf(G, J-K, N-J), K-2, "\n", sep=", ")
  
  p3.T <- update(p2.T, . ~ . + I(emcmids^3))
  clustersub <- emcp[attributes(p3.T$model)$row.names]
  c3<-clx(p3.T, 1,clustersub)
  K <- p3.T$rank 
  ESSr <- sum(p3.T$residuals^2)
  G = ((ESSr-ESSur)/(J-K)) / (ESSur / (N-J)) # ~F(J-K, N-J)
  cat(pval <- pf(G, J-K, N-J), K-2, "\n", sep=", ")

  p4.T <- update(p3.T, . ~ . + I(emcmids^4))
  clustersub <- emcp[attributes(p4.T$model)$row.names]
  clx(p4.T, 1, clustersub)
  K <- p4.T$rank 
  ESSr <- sum(p4.T$residuals^2)
  G = ((ESSr-ESSur)/(J-K)) / (ESSur / (N-J)) # ~F(J-K, N-J)
  cat(pval <- pf(G, J-K, N-J), K-2, "\n", sep=", ")
  
  p5.T <- update(p4.T, . ~ . + I(emcmids^5))
  clustersub <- emcp[attributes(p5.T$model)$row.names]
  clx(p5.T, 1, clustersub)
  K <- p5.T$rank 
  ESSr <- sum(p5.T$residuals^2)
  G = ((ESSr-ESSur)/(J-K)) / (ESSur / (N-J)) # ~F(J-K, N-J)
  cat(pval <- pf(G, J-K, N-J), K-2, "\n", sep=", ")  
  
  p6.T <- update(p5.T, . ~ . + I(emcmids^6))
  clustersub <- emcp[attributes(p6.T$model)$row.names]
  clx(p6.T, 1, clustersub)
  K <- p6.T$rank 
  ESSr <- sum(p6.T$residuals^2)
  G = ((ESSr-ESSur)/(J-K)) / (ESSur / (N-J)) # ~F(J-K, N-J)
  cat(pval <- pf(G, J-K, N-J), K-2, "\n", sep=", ")  
  
  p7.T <- update(p6.T, . ~ . + I(emcmids^7))
  clustersub <- emcp[attributes(p7.T$model)$row.names]
  clx(p7.T, 1, clustersub)
  K <- p7.T$rank 
  ESSr <- sum(p7.T$residuals^2)
  G = ((ESSr-ESSur)/(J-K)) / (ESSur / (N-J)) # ~F(J-K, N-J)
  cat(pval <- pf(G, J-K, N-J), K-2, "\n", sep=", ")
 
  ################################### Global fits end

  ## regs run on collapsed to cell level data with heteroskedasticity-consistent standard errors.
  ##Also include variance correction factor as in Card and Lee 2008 equations 12 and 13.
  ## Omit for now
  
  #p0.C <- lm(binned.outcome~forcing, weights=binned.takers)
  #sqrt(vcovHC(p0.C, "HC")) #robust Heteroskedasticity-consistent standard errors: about equal to clx(p0.T, 1, clustersub)
  #sigmahat2_a <- 2*(sum(p0.C$residuals^2)/sum(Nj) - sum(binned.variance, na.rm=T)/sum(Nj))  #maybe this one?
  #newstderror <- vcovHC(p0.C, "HC")[2,2] + sigmahat2_a

  #p1.C <- lm(binned.outcome~forcing*breaks[2:81], weights=binned.takers)
  #sqrt(vcovHC(p1.C, "HC"))
  #sigmahat2_a <- 2*(sum(p1.C$residuals^2 * Nj[1:79]) / sum(Nj)  - sum(binned.variance, na.rm=T)/sum(Nj)) #maybe this one?
  #newstderror <- vcovHC(p1.C, "HC")[2,2] + sigmahat2_a

  #p2.C <- lm(outcome~forcing.sample*emc + forcing.sample*I(emc^2), subset=obs==1)
  return(list(c1, c2, c3))
}


srdregs <- function(binwidth, maxbin, outcome, testtype, emcp=emc, coefs=1, use.coefs=FALSE) {
    ## outcome = const + B1*age + B2*Cutoff + e
    ## outcome = const + B1*emcp + B2*forcing.sample + e

  breaks <- getbins(binwidth, maxbin)@breaks
  emcdum <- ifelse(emcp >= min(breaks) & emcp <= max(breaks), 1, 0)

  cutbybreak <- cut(emcp[emcdum==1], breaks=breaks)
  binmids <- hist(emcp[emcdum==1], breaks=breaks, plot=FALSE)$mids

  emcbinmids <- findInterval(emcp[emcdum==1], binmids, rightmost.closed=TRUE)
  emcmids <-  binmids[emcbinmids]
  
  forcing <- ifelse(binmids>0, 1, 0)
  #tx <- which(forcing==1)
  forcing.sample <- ifelse(emcp>0, 1, 0)

  test <- outcome[emcdum==1][emcbinmids!=0]
  alpha <- forcing.sample[emcdum==1][emcbinmids!=0]
  if (use.coefs){
    beta <- coefs[emcdum==1,][emcbinmids!=0,]
    sr.model <- lm(test~alpha+emcmids+I(emcmids^2)+beta)
  } else sr.model <- lm(test~alpha+emcmids+I(emcmids^2))
  
  clustersub <- emcp[attributes(sr.model$model)$row.names]
  adjrsq <- summary(sr.model)$adj.r.squared

  ##temporary fix for now until I figure out prob with sandwich estimator
  a<-try(print(clx(sr.model, 1, clustersub)), silent=TRUE)
  if (class(a)=="try-error") print(summary(sr.model))
  cat("N=",length(clustersub), " Adj. R-sq=", adjrsq, '\n\n', sep='')
  

}

frdregs <- function(binwidth, maxbin, outcome, testtype, emcp=emc, coefs=1, use.coefs=FALSE, instr) {
    ## outcome = const + B1*age + B2*Cutoff + e
    ## outcome = const + B1*emcp + B2*forcing.sample + e
    ## forcing.sample =  const+A1*age

  breaks <- getbins(binwidth, maxbin)@breaks
  emcdum <- ifelse(emcp >= min(breaks) & emcp <= max(breaks), 1, 0)

  cutbybreak <- cut(emcp[emcdum==1], breaks=breaks)
  binmids <- hist(emcp[emcdum==1], breaks=breaks, plot=FALSE)$mids

  emcbinmids <- findInterval(emcp[emcdum==1], binmids, rightmost.closed=TRUE)
  emcmids <-  binmids[emcbinmids]
  
  forcing <- ifelse(binmids>0, 1, 0)
  #tx <- which(forcing==1)
  forcing.sample <- ifelse(emcp>0, 1, 0)

  test <- outcome[emcdum==1][emcbinmids!=0]
  alpha <- forcing.sample[emcdum==1][emcbinmids!=0]
  instrument <- instr[emcdum==1][emcbinmids!=0]
  
  if (use.coefs){
    beta <- coefs[emcdum==1,][emcbinmids!=0,]
    mf <- data.frame(test, alpha, instrument, emcmids, beta)
    form1 <- as.formula(paste("test~emcmids+I(emcmids^2)+alpha+",
                   paste(colnames(mf[,5:dim(mf)[2]]), collapse="+"),
                         " | instrument+emcmids+I(emcmids^2)+",paste(colnames(mf[,5:dim(mf)[2]]), collapse="+"), sep=''))
    sr.model <- ivreg(form1, data=mf)
  } else {
    mf <- data.frame(test, alpha, instrument, emcmids)
    sr.model <- ivreg(test~emcmids+I(emcmids^2)+alpha | instrument+emcmids+I(emcmids^2), data=mf)
  }

  clustersub <- emcp[attributes(sr.model$model)$row.names]
  adjrsq <- summary(sr.model)$adj.r.squared

  ##temporary fix for now until I figure out prob with sandwich estimator
  a<-try(print(clx(sr.model, 1, clustersub)), silent=TRUE)
  if (class(a)=="try-error") print(summary(sr.model))
  cat("N=",length(clustersub), " Adj. R-sq=", adjrsq, '\n\n', sep='')
  

}
 
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ FUNCTIONS END $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$##
