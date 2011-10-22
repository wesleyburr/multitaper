##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2011 Karim Rahim 

##     This file is part of the multitaper package for R.

##     The multitaper package is free software: you can redistribute it and
##     or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, either version 2 of the License, or
##     any later version.

##     The multitaper package is distributed in the hope that it will be 
##     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
##     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.

##     You should have received a copy of the GNU General Public License
##     along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

##     If you wish to report bugs please contact the author. 
##     karim.rahim@gmail.com
##     112 Jeffery Hall, Queen's University, Kingston Ontario
##     Canada, K7L 3N6



## these are plotting utilities. Usage for other purposes is not supported.

## dropFreqs Class for selecting frequencies for plotting
## can be used for spec, mtm and mtm.coh classes.

dropFreqs <- function(spec, minFreq, maxFreq) UseMethod("dropFreqs")

dropFreqs.default <- function(spec, minFreq, maxFreq) {
    print("This function is only valid for objects of spec, mtm, or mtm.coh classes")
    spec
}

dropFreqs.spec <- function(spec, minFreq, maxFreq) {
    idx <- (findInterval(spec$freq, c(minFreq,maxFreq)) == 1)

    if(sum(idx) <= 1) {
        stop("minFreq and maxFreq must allow for a range of frequencies to be returned")
    }
    spec.out <- spec
    spec.out$freq <- spec$freq[idx]
    spec.out$spec <- spec$spec[idx]

    spec.out
}

dropFreqs.mtm <- function(spec, minFreq, maxFreq) {
    idx <- (findInterval(spec$freq, c(minFreq,maxFreq)) == 1)

    if(sum(idx) <= 1) {
        stop("minFreq and maxFreq must allow for a range of frequencies to be returned")
    }
    spec.out <- spec
    spec.out$freq <- spec$freq[idx]
    spec.out$spec <- spec$spec[idx]

    ##adjust mtm parameters
    if(!is.null(spec.out$mtm)) {
        ## null unnecessary values 
        ## enforces fact that currently function is mainly a
        ## plotting utility
        spec.out$mtm$dpss <- NULL
        spec.out$mtm$eigenCoefs <- NULL
        spec.out$mtm$eigenCoefsWt <- NULL
        
        ## keep values used in plotting
        spec.out$mtm$Ftest <- spec.out$mtm$Ftest[idx]
        spec.out$mtm$dofs <- spec.out$mtm$dofs[idx]

        if(!is.null(spec.out$mtm$jk)) {
            spec.out$mtm$jk$varjk <- NULL
            spec.out$mtm$jk$upperCI <- spec.out$mtm$jk$upperCI[idx]
            spec.out$mtm$jk$maxVal <- max(spec.out$mtm$jk$upperCI)
            spec.out$mtm$jk$bcjk <- NULL
            spec.out$mtm$jk$lowerCI <- spec.out$mtm$jk$lowerCI[idx]
            spec.out$mtm$jk$sjk <- NULL
            spec.out$mtm$jk$minVal <- min(spec.out$mtm$jk$lowerCI)
        }
    }

    spec.out
}

dropFreqs.mtm.coh <- function(spec, minFreq, maxFreq) {
    idx <- (findInterval(spec$freq, c(minFreq,maxFreq)) == 1)

    if(sum(idx) <= 1) {
        stop("minFreq and maxFreq must allow for a range of frequencies to be returned")
    }

    spec.out <- spec
    spec.out$NTmsc <- spec.out$NTmsc[idx]
    spec.out$msc <- spec.out$msc[idx]
    spec.out$NTvar <- spec.out$NTvar[idx]
    spec.out$freq <- spec.out$freq[idx]
    spec.out$ph <- spec.out$ph[idx]
    spec.out$phvar <- spec.out$phvar[idx]
    spec.out$nfreqs <- sum(idx)

    spec.out
}
