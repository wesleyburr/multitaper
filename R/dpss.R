##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2011 Karim Rahim 
##
##     Written by Karim Rahim, with code from David J. Thomson.
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

##########################################################
##
##  dpss
## 
##  Generates k orthogonal discrete prolate spheroidal
##  sequences (dpss) using the tridiagonal method. See
##  Thomson, D.J., Spectrum Estimation and Harmonic Analysis,
##  Proceedings of the IEEE, 1982. 
##
##########################################################

dpss <- function(n, k, nw, returnEigenvalues=TRUE) {

    stopifnot(n >= 1, nw >= 0.5, k >= 1, nw <= 500, k <= 1.5+2*nw) 
        
    ##eigen is of length k for use by LAPACK functions.
    out <- .Fortran("dpss", as.integer(n), as.integer(k),
              as.double(nw), 
              v=double(n*k), eigen=double(k),
              PACKAGE='multitaper')

    out$v <- matrix(out$v, nrow=n, ncol=k, byrow=FALSE)
    if(returnEigenvalues) {
        out$eigen <- dpssToEigenvalues(out$v, nw)
    } else {
        out$eigen <- NULL
    }
    
    res <- list(v=out$v,
                eigen=out$eigen)
    class(res) <- "dpss"
    return(res)
}

##########################################################
##
##  dpssToEigenvalues
## 
##  Given a set of dpss tapers, find the eigenvalues after
##  the fact. 
##
##########################################################

dpssToEigenvalues <- function(v, nw) {
    v <- as.matrix(v)
    n <- length(v[,1])
    k <- length(v[1,])
    w <- nw/n
    npot <- 2**(ceiling(log2(2*n)))

    ## pad
    scratch <- rbind(v, matrix(data=0, nrow=npot-n, ncol=k))
    ## n * acvs
    scratch <- Re(mvfft(abs(mvfft(scratch))**2))/npot
    
    j <- 1:(n-1)
    vectorOfRatios <- c(2*w, sin(2*pi*w*j)/(pi*j))
    
    ## Note: both vector of ratios and scratch
    ##  roughy decrease in amplitude with increasing index,
    ##  so sum things up in reverse order
    eigenvalues <- NULL
    if(k>1) {
        eigenvalues <- apply(scratch[n:2,] * vectorOfRatios[n:2], 2, sum)
    } else {
        eigenvalues <- sum(scratch[n:2,] * vectorOfRatios[n:2])
    }
    return(2*eigenvalues +
           vectorOfRatios[1]*scratch[1,1:k])
}

