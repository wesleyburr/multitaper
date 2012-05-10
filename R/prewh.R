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


##############################################################
##
##  prewh
##
##  Prewhiten data using Spectrum->ACV as prewhitener.
##
##############################################################
prewh <- function(x, dT, maxlag){

  mX <- mean(x)
  x <- x - mX
  Sp <- spec.mtm(x,dT=dT,plot=FALSE)
  acv <- Sp2Acv(Sp,maxlag)

  N <- length(x)
  y <- rep(NA,N)
  coef <- -1*acv[(maxlag+1):1]/acv[1]
  coef[(maxlag+1)] <- abs(coef[(maxlag+1)])

  for(j in (maxlag+1):N) {
    y[j] <- sum(coef * x[(j-maxlag):j])
  }
  y 
}

##############################################################
##
##  Sp2Acv
##
##  Takes spectrum estimate ('spec' object) and inverts it
##  to find a robust (if mtm is used) estimate of the ACV.
##
##############################################################
Sp2Acv <- function(Sp,maxlag) {

    # requires a spec object; works with mtm objects as well
  if(!("spec" %in% class(Sp))) {
    stop("Sp2Acv requires a 'spec' object to function.")
  }

  if("mtm" %in% class(Sp)) {
    nFFT <- Sp$mtm$nFFT  
    nFrHi <- Sp$mtm$nfreqs
    x <- rep(0,nFFT)
    y <- x
    if((nFFT %% 2 == 0) & (nFrHi %% 2 != 0)) {
      x[1:nFrHi] <- Sp$spec
      x[nFFT:(nFrHi+1)] <- Sp$spec[2:(nFrHi-1)]
    } else {
      nFrHi <- nFrHi + 1
      x[1] <- 0.0
      x[2:nFrHi] <- Sp$spec
      x[nFFT:(nFrHi+1)] <- Sp$spec[1:(nFrHi-2)]
    }
    delF <- Sp$freq[2] - Sp$freq[1]
  } else {
    # default spectrum object does not contain zero frequency
    nFFT <- length(Sp$freq)*2
    nFrHi <- length(Sp$freq)+1
    delF <- Sp$freq[2] - Sp$freq[1]
    x <- rep(0,nFFT)
    y <- x
    x[1] <- 0.0
    x[2:nFrHi] <- Sp$spec
    x[nFFT:(nFrHi+1)] <- Sp$spec[1:(nFrHi-2)]
  }
  x <- x*delF
  acv8 <- Re(fft(complex(real=x,imaginary=y),inverse=TRUE))
  acv8[1:(maxlag+1)]
}

