##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2011 Karim Rahim 
##
##     Written by Karim Rahim, with contributions from Wesley Burr.
##
##     This file is part of the multitaper package for R.
##     http://cran.r-project.org/web/packages/multitaper/index.html
## 
##     The multitaper package is free software: you can redistribute it and 
##     or modify it under the terms of the GNU General Public License as 
##     published by the Free Software Foundation, either version 2 of the 
##     License, or any later version.
##
##     The multitaper package is distributed in the hope that it will be 
##     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
##     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License
##     along with multitaper.  If not, see <http://www.gnu.org/licenses/>.
##
##     If you wish to report bugs please contact the author:
## 
##     Karim Rahim
##     karim.rahim@gmail.com
##     112 Jeffery Hall, Queen's University, Kingston Ontario
##     Canada, K7L 3N6

spec.mtm <- function(timeSeries,
                     nw=5.0,
                     k=10,
                     nFFT="default", 
                     taper=c("dpss"),
                     centre=c("Slepian"),
                     dpssIN=NULL,
                     returnZeroFreq=TRUE,
                     Ftest=FALSE,
                     jackknife=FALSE,
                     jkCIProb=.95,
                     maxAdaptiveIterations=100,
                     plot=TRUE,
                     na.action=na.fail,
                     returnInternals=FALSE,
                     ...) {

    series <- deparse(substitute(timeSeries))
    taper <- match.arg(taper,c("dpss","sine"))
    centre <- match.arg(centre,c("Slepian","arithMean","trimMean","none"))

    if( (taper=="sine") && jackknife) { stop("Cannot jackknife over sine tapers.")}
    if( (taper=="sine") && Ftest) { stop("Cannot compute Ftest over sine tapers.")}

    # warning for deltaT missing: makes all frequency plots incorrect
    if(!is.ts(timeSeries)) {
      warning("Time series is not a ts object. deltaT is not set, and frequency axes may be incorrect.")
      deltaT <- 1
      timeSeries <- as.double(as.ts(timeSeries))
    } else {
      deltaT <- deltat(timeSeries)
      timeSeries <- as.double(timeSeries)
    }
    n <- length(timeSeries)

    if(taper=="dpss") {
      stopifnot(nw >= 0.5, k >= 1, nw <= 500, k <= 1.5+2*nw, n > 8)
    } else {
      stopifnot(k <= n, k >= 1, n > 8)
    }

    na.action(timeSeries)
    if(k==1) {
        Ftest=FALSE
        jackknife=FALSE
    }

    sigma2 <- var(timeSeries) * (n-1)/n

    if(nFFT == "default") {
        nFFT <- 2* 2^ceiling(log2(n))
    } else {
        stopifnot(is.numeric(nFFT))
    }
    stopifnot(nFFT >= n)

    ## convert time-series to zero-mean by one of three methods, if set; default is Slepian
    if(centre=="Slepian") {
      timeSeries <- centre(timeSeries, nw=nw, k=k, deltaT=deltaT)    
    } else if(centre=="arithMean") {
      timeSeries <- centre(timeSeries, trim=0) 
    } else if(centre=="trimMean") {
      timeSeries <- centre(timeSeries, trim=0.10)
    }

    if(taper=="dpss") { 
      mtm.obj <- spec.mtm.dpss(timeSeries=timeSeries,
                     nw=nw, k=k, nFFT=nFFT, 
                     dpssIN=dpssIN, returnZeroFreq=returnZeroFreq, 
                     Ftest=Ftest, jackknife=jackknife, jkCIProb=jkCIProb, 
                     maxAdaptiveIterations=maxAdaptiveIterations, 
                     returnInternals=returnInternals, 
                     n=n, deltaT=deltaT, sigma2=sigma2, series=series,
                     ...) 
    } else if(taper=="sine") {
      mtm.obj <- spec.mtm.sine(timeSeries=timeSeries,
                     nFFT=nFFT, dpssIN=dpssIN, returnZeroFreq=returnZeroFreq,
                     returnInternals=FALSE, n=n, deltaT=deltaT, sigma2=sigma2,
                     series=series, ...)
    }

    if(plot) {
        plot.mtm(mtm.obj, jackknife=jackknife, ...)
        return(invisible(mtm.obj))
    } else {
        return(mtm.obj)
    }
}

##############################################################
##
##  spec.mtm.dpss
##
##  Computes multitaper spectrum using Slepian tapers
##
##############################################################
spec.mtm.dpss <- function(timeSeries,
                     nw,
                     k,
                     nFFT,
                     dpssIN,
                     returnZeroFreq,
                     Ftest,
                     jackknife,
                     jkCIProb,
                     maxAdaptiveIterations,
                     returnInternals,
                     n,
                     deltaT,
                     sigma2,
                     series,
                     ...) {

    dw <- NULL
    ev <- NULL
    receivedDW <- TRUE
    if(!.is.dpss(dpssIN)) {
      receivedDW <- FALSE
      dpssIN <- dpss(n, k, nw=nw, returnEigenvalues=TRUE)
      dw <- dpssIN$v*sqrt(deltaT)
      ev <- dpssIN$eigen 
    }
    else {
      dw <- .dpssV(dpss)
      ev <- .dpssEigen(dpss)
      if(ev == NULL) {
        ev <- dpssToEigenvalues(dw, nw) }
        dw <- dw*sqrt(deltaT) 
    }

    nFreqs <- nFFT %/% 2 + as.numeric(returnZeroFreq)
    offSet <- if(returnZeroFreq) 0 else 1 
    scaleFreq <- 1 / as.double(nFFT * deltaT)
    
    swz <- NULL ## Percival and Walden H0
    ssqswz <- NULL

    taperedData <- dw*timeSeries
    
    nPadLen <- nFFT - n
    paddedTaperedData <- rbind(taperedData, matrix(0, nPadLen, k))
    cft <- mvfft(paddedTaperedData)
    cft <- cft[(1+offSet):(nFreqs+offSet),]
    sa <- abs(cft)**2
    
    resultFreqs <- ((0+offSet):(nFreqs+offSet-1))*scaleFreq 

    adaptive <-  NULL
    jk <- NULL
    PWdofs <- NULL
    if(!jackknife) {
        adaptive <- .mw2wta(sa, nFreqs, k,
                            sigma2, deltaT, ev)
    } else {
        stopifnot(jkCIProb < 1, jkCIProb > .5)
        adaptive <- .mw2jkw(sa, nFreqs, k,
                            sigma2, deltaT, ev)
        scl <- exp(qt(jkCIProb,adaptive$dofs)*
                   sqrt(adaptive$varjk))
        upperCI <- adaptive$s*scl
        lowerCI <- adaptive$s/scl
        minVal = min(lowerCI)
        maxVal = max(upperCI)
        jk <- list(varjk=adaptive$varjk,
                   bcjk=adaptive$bcjk,
                   sjk=adaptive$sjk,
                   upperCI=upperCI,
                   lowerCI=lowerCI,
                   maxVal=maxVal,
                   minVal=minVal)
    }

    ftestRes <- NULL

    if(Ftest) {
        if(is.null(swz)) {
            swz <- apply(dw, 2, sum)
        }
        ftestRes <- .HF4mp1(cft,
                            swz,
                            k,
                            ssqswz)
    }

    eigenCoef1 <- NULL
    wtCoef1 <- NULL
    
    if(returnInternals) {
        eigenCoef1 <- cft
        wtCoef1 <- sqrt(adaptive$wt)
    }
    auxiliary <- list(dpss=dpssIN,
                      eigenCoefs=eigenCoef1,
                      eigenCoefWt=wtCoef1,
                      nfreqs=nFreqs,
                      nFFT=nFFT,
                      jk=jk,
                      Ftest=ftestRes$F,
                      dofs=adaptive$dofs,
                      nw=nw,
                      k=k)

    spec.out <- list(origin.n=n,
                     method="Multitaper Spectral Estimate",
                     pad= nFFT - n,
                     spec=adaptive$s,
                     freq=resultFreqs,
                     series=series,
                     mtm=auxiliary)

    class(spec.out) <- c("mtm", "spec")
    
    if(Ftest) {
        class(spec.out) <- c("mtm", "Ftest", "spec")
    }
    return(spec.out);
}


#########################################################################
##
##  spec.mtm.sine 
##
##  Computes multitaper spectrum estimate using sine tapers, as in
## 
##  Riedel, Kurt S. and Sidorenko, Alexander, Minimum Bias Multiple 
##    Taper Spectral Estimation. IEEE Transactions on Signal Processing,
##    Vol. 43, No. 1, January 1995.
## 
#########################################################################

spec.mtm.sine <- function(timeSeries,
                          nFFT,
                          dpssIN, 
                          returnZeroFreq,
                          n, 
                          deltaT, 
                          sigma2,
                          series=series,
                          ...) {
                          
    dw <- NULL
    receivedDW <- TRUE
    if(!.is.dpss(dpssIN)) {
      receivedDW <- FALSE
      dpssIN <- sineTaper(n, k)
      dw <- dpssIN$v
    }
    else {
      dw <- .dpssV(dpss)
    }

    nFreqs <- nFFT %/% 2 + as.numeric(returnZeroFreq)
    offSet <- if(returnZeroFreq) 0 else 1 
    scaleFreq <- 1 / as.double(nFFT * deltaT)
    
    taperedData <- dw*timeSeries
    
    nPadLen <- nFFT - n
    paddedTaperedData <- rbind(taperedData, matrix(0, nPadLen, k))
    cft <- mvfft(paddedTaperedData)
    cft <- cft[(1+offSet):(nFreqs+offSet),]
    
    resultFreqs <- ((0+offSet):(nFreqs+offSet-1))*scaleFreq 

    ## compute sine-based multitaper estimate
    ## cft is nFreqs by k
    ## Formula from Abstract, Riedel & Sidorenko

    ## first find y(f \pm j/(2n+2) ) for j=1,2,...,k
    ##  --> requires either wrapping around 0/Nyquist, or ignoring missing elements

    ## set up return object
    auxiliary <- list(dpss=dpssIN,
                      eigenCoefs=NULL,
                      eigenCoefWt=NULL,
                      nfreqs=nFreqs,
                      nFFT=nFFT,
                      jk=NULL,
                      Ftest=NULL,
                      dofs=NULL,
                      nw=NULL,
                      k=k)

    spec.out <- list(origin.n=n,
                     method="Sine Multitaper Spectral Estimate",
                     pad= nFFT - n,
                     spec=Sfinal,
                     freq=resultFreqs,
                     series=series,
                     mtm=auxiliary)

    class(spec.out) <- c("mtm", "spec")
    return(spec.out);
}


#########################################################################
##
## centre
##
## Takes a time series and converts to zero-mean using one of three 
## methods: Slepian projection, arithmetic mean, or trimmed mean.
## 
#########################################################################

centre <- function(x, nw=NULL, k=NULL, deltaT=NULL, trim=0) {
    na.fail(x)
    res <- NULL
    if(is.null(nw) && is.null(k) ) {
        res <- x - mean(x, trim=trim)
    } else {
        if(trim != 0) {
            warning(paste("Ignoring trim =", trim))
        }
        stopifnot(nw >= 0.5, k >= 1, nw <= 500, k <= 1.5+2*nw)
        if(is.null(deltaT)) {
            if(is.ts(x)) {
                deltaT <- deltat(ts)
            } else {
                warning("deltaT not specified; using deltaT=1.")
                deltaT <- 1
            }
        }
        n <- length(x)
        dpssRes <- dpss(n, k=k, nw=nw,
                        returnEigenvalues=TRUE)
        dw <- dpssRes$v*sqrt(deltaT)
        ev <- dpssRes$eigen
        swz <- apply(dw, 2, sum)
        ## zero swz where theoretically zero
        swz[1:(k/2)*2] <- 0.0
        ssqswz <- sum(swz**2)
        res <- .mweave(x, dw, swz,
                       n, k, ssqswz, deltaT)
        res <- x - res$cntr
    }
    return(res)
}


#########################################################################
##
## jackknife coherence and helper smoother and plotting functions
## 
## Example: 
## jkRes <- jkcoh1(r1$auxiliary$cft, r2$auxiliary$cft,
##                 4,2048,4,4096,395)
## pGreater <-  percentjkMSCGreaterThan(jkRes$msc, 4)
## plotJkcoh1(r1$freqs, jkRes$TRmsc, jkRes$NTvar, 4, pGreater)
##
#########################################################################

mtm.coh <- function(mtm1, mtm2,
                   fr=NULL, tau=0, plot=TRUE,...) {

    ## note Dave saves the cft
    ## in ./odinlibs-1.1/src/mw/mw2pakt as weighted
    ## 1000 blkcft(n,k,curblk,curset) =
    ##  cft(n*ndecfr,k)*sqrt(wt(n*ndecfr,k))

    ## we require auxiliary data
    if(is.null(mtm1$mtm$eigenCoefs) || is.null(mtm2$mtm$eigenCoefs)) {
        stop("Both mtm objects must have been computed with returnInternals=TRUE.")
    }

    if(mtm1$mtm$k != mtm1$mtm$k) {
        stop("Both mtm objects must have the same value for k.")
    }
    ##k <- mtm1$auxiliary$

    if(mtm1$mtm$nfreqs != mtm1$mtm$nfreqs) {
        stop("Both mtm objects must have the same value for nFFT.")
    }

    nord <- mtm1$mtm$k
    nfreqs <- mtm1$mtm$nfreqs
    cft1 <- mtm1$mtm$eigenCoefs
    cft2 <- mtm2$mtm$eigenCoefs
    
    fr <-  if(is.null(fr))  array(as.double(0), nfreqs) else fr
    
    blklof <-  if(nfreqs %%2 ==0) 1 else 0
    blkhif <- nfreqs -1 + blklof

    nordP2 <- nord +2
    phcorr <-  1
    out <- .Fortran("jkcoh1", cft1=as.complex(cft1),
                    cft2=as.complex(cft2), nord=as.integer(nord),
                    blklof=as.integer(blklof), blkhif=as.integer(blkhif),
                    fr=as.double(fr),  tau=as.double(tau),
                    phcorr=as.integer(phcorr),
                    NTmsc=double(nfreqs), NTvar=double(nfreqs),
                    msc=double(nfreqs), ph=double(nfreqs),
                    phvar=double(nfreqs),
                    s1=double(nordP2), s2=double(nordP2),
                    jkmsc=double(nordP2), TRmsc=double(nordP2),
                    bias=double(nfreqs),
                    cx=complex(nordP2))

    coh.out <- list(NTmsc=out$NTmsc, NTvar=out$NTvar,
                    msc=out$msc, nfreqs=mtm1$mtm$nfreqs,
                    freq=mtm1$freq, k=nord,
                    ph=out$ph, phvar=out$phvar)
    class(coh.out) <- "mtm.coh"
    
    
   if(plot) {
        plot.mtm.coh(coh.out, ...)
        return(invisible(coh.out))
    } else {
        return(coh.out)
    }
}


