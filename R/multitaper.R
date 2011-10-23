##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2011 Karim Rahim 
##
##     Written by Karim Rahim and Wesley Burr.
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
                     sineAdaptive=FALSE,
                     ...) {

    series <- deparse(substitute(timeSeries))
    taper <- match.arg(taper,c("dpss","sine"))
    centre <- match.arg(centre,c("Slepian","arithMean","trimMean","none"))

    if( (taper=="sine") && jackknife) { stop("Cannot jackknife over sine tapers.")}
    if( (taper=="sine") && Ftest) { stop("Cannot compute Ftest over sine tapers.")}
    if( (taper=="sine") && !returnZeroFreq) { returnZeroFreq = TRUE; 
                                              warning("returnZeroFreq must be TRUE for sine taper option.") }

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

    # *** should clean this up; we compute swz and ssqswz in centre(), then have to recompute
    #  it below for the Ftest ...
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
      mtm.obj <- spec.mtm.sine(timeSeries=timeSeries, k=k, sineAdaptive=sineAdaptive,
                     nFFT=nFFT, dpssIN=dpssIN, returnZeroFreq=returnZeroFreq,
                     returnInternals=FALSE, n=n, deltaT=deltaT, sigma2=sigma2,
                     series=series,maxAdaptiveIterations=maxAdaptiveIterations,
                     ...)
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
    swz <- apply(dw, 2, sum)
    swz[1:(k/2)*2] <- 0.0
    ssqswz <- sum(swz**2)

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
##  Algorithm implementation based on previous work by:
##    German Prieto, Universidad de los Andes
##       via \texttt{mtsepc}, a F90 package that can be found at
##       http://wwwprof.uniandes.edu.co/~gprieto/software/mwlib.html
##
##    and
##
##    Robert L. Parker, Scripps Institution of Oceanography
##      via \texttt{psd.f}, a F77 program that can be found at
##      http://igppweb.ucsd.edu/~parker/Software/Source/psd.f
## 
#########################################################################

spec.mtm.sine <- function(timeSeries,
                          nFFT,
                          k,
                          sineAdaptive,
                          dpssIN, 
                          returnZeroFreq=TRUE,
                          n, 
                          deltaT, 
                          sigma2,
                          series=series,
                          maxAdaptiveIterations,
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

    nFFT <- nFFT*2
    # returnZeroFreq forced to TRUE, offset = 0
    nFreqs <- nFFT %/% 4 + as.numeric(returnZeroFreq)
    offSet <- if(returnZeroFreq) 0 else 1 
    scaleFreq <- 1 / as.double(nFFT/2 * deltaT)
    resultFreqs <- ((0+offSet):(nFreqs+offSet-1))*scaleFreq 
    nPadLen <- nFFT - n
    df <- 1/nFFT/deltaT

    # compute a single FFT; since we are using sine tapers, this is all we need
    ones <- matrix(1,n,1)
    timeSeries <- timeSeries*ones
    paddedData<- rbind(timeSeries, matrix(0, nPadLen, 1))
    cft <- mvfft(paddedData)
 
    # constant number of tapers, or adaptive?
    spec <- as.double(matrix(0,1,nFreqs))

    if(!sineAdaptive) { # constant k tapers
      spec <- quickSine(nFreqs=nFreqs,nFFT=nFFT,k=k,cft=cft,useAdapt=FALSE,kadapt=NULL)      
    } else {
      # smoothing factor defaults to 1, not changeable
      fact <- 1;
      initTaper <- ceiling(3.0 + sqrt(fact*n)/5.0);
 
      # c_1, c_2 are constants for parabolic weighting
      c1=1.2000 
      c2=3.437

      # pilot estimate of S
      spec0 <- quickSine(nFreqs=nFreqs,nFFT=nFFT,k=initTaper,
                         cft=cft,useAdapt=FALSE,kadapt=NULL)

      # initialize kadapt
      kadapt <- matrix(data=k, nrow=nFreqs, ncol=1)
      opt <- matrix(0.0, nrow=nFreqs, ncol=1)
      for(j in 1:maxAdaptiveIterations) {
        #cat(".")
        y <- log(spec0);

        system.time(
        for(m in 1:nFreqs) {
            # Estimate kadapt(f) for each f
            # R = S''/S ~ y'' + (y')^2 for y=ln(S)
  
            ispan = round(kadapt[m]*1.4)
            deriv <- northog(n=nFreqs,i1=(m-ispan),i2=(m+ispan),s=y)
            dy <- deriv[1]
            ddy <- deriv[2]
  
            R <- (ddy  + dy^2)/df^2
            ak <- kadapt[m]/(2*ispan) # correct for integer steps 
            phi <- 720.0 * ak^5 * (1.0 - 1.286*ak + 0.476*ak^3 - 0.0909*ak^5)
            sigR <- sqrt(phi/(kadapt[m]^5)) / df^2
            opt[m] <- c2/(df^4 *( R^2 + 1.4*sigR^2) /fact^2)^0.2
          } # end of frequencies
        )
        opt <- curb(n=nFreqs,v=opt)
        opt[opt <= 3] <- 3
        kadapt <- opt
        # recompute spectra for next step in iteration
        spec0 <- quickSine(nFreqs=nFreqs,nFFT=nFFT,cft=cft,useAdapt=TRUE,kadapt=kadapt)
      } # end of iterative loop    
      #cat("\n");
      spec <- spec0;
    } # end of adaptive logic

    # normalize spectrum
    const <- var(timeSeries)/sum(spec)/df
    specFinal <- const*spec

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
                     spec=specFinal,
                     freq=resultFreqs,
                     series=series,
                     mtm=auxiliary)

    class(spec.out) <- c("mtm", "spec")
    return(spec.out);
}

#########################################################################
##
## quickSine
##
## Sine Taper constant taper number quick iterative spectrum estimation
## Based on same source as spec.mtm.sine.
## 
#########################################################################

quickSine <- function(nFreqs,nFFT,k,cft,useAdapt,kadapt) { 

      spec <- rep(0.0, nFreqs)
      for(i in 1:nFreqs) {
        i2 <- 2*(i-1);
        if(useAdapt) {
          ks <- kadapt[i] 
        } else {
          ks <- k
        }
        ck <- 1/(ks^2);

        # Use parabolic weighting
        for(j in 1:ks) {
          j1 <- (i2+nFFT-j) %% nFFT;
          j2 <- (i2+j) %% nFFT;
          zz <- cft[j1+1] - cft[j2+1];
          wt <- 1.0 - ck*(j-1)^2;

          spec[i] = spec[i] + (Mod(zz)^2) * wt;
        }
        # normalize for parabolic factor
        spec[i] = spec[i] * (6.0*ks)/(4*(ks^2)+(3*ks)-1);
      } 
      return(spec)
} 


#########################################################################
##
## northog
##
## Performs quadratically-weighted LS fit to some function 's' by
## a degree-two polynomial in an orthogonal basis; returns
## d1 and d2, estimates of 1st and 2nd derivatives at center of record
## 
#########################################################################


northog <- function(n,i1,i2,s) {
      L = i2 - i1 + 1
      el=L
      gamma = (el^2 - 1.0)/12.0
      u0sq = el
      u1sq = el*(el^2 - 1.0)/12.0
      u2sq = (el*(el^2 - 1.0)*(el^2- 4.0))/180.0
      amid= 0.5*(el + 1.0)
      dot0=0.0
      dot1=0.0
      dot2=0.0
      ssq=0.0
      for(kk in 1:L) {
        i=kk + i1 - 1
        if (i <= 0) { i=2 - i;}
        if (i > n) { i=2*n - i;}
        dot0 = dot0 + s(i)
        dot1 = dot1 + (kk - amid) * s(i)
        dot2 = dot2 + ((kk - amid)^2 - gamma)*s(i)
      }
      ds = dot1/u1sq
      dds = 2.0*dot2/u2sq
      return(c(ds,dds))
}

#########################################################################
##
## curb
##
##  Reworks the input n-vector v() so that all points lie below
##  the piece-wise linear function v(k) + abs(j-k), where v(k)
##  is a local minimum in the original v.
##  Effectively clips strong peaks and keeps slopes under 1 in
##  magnitude.
## 
#########################################################################
curb <- function(n, vin) {
      v <- vin;
      for(j in 2:(n-1)) {
        if (v[j] < v[j+1] && v[j] < v[j-1]) { vloc <- v[j]; 
          for(k in 1:n) {
            v[k] <- min(v[k], vloc+abs(j-k));
          }
        }
      }
      return(v[1:n]);
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


