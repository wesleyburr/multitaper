##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2011 Karim Rahim 
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

.mweave <-  function (x,dw,swz,ndata,nord,ssqswz,dt_) {
    
    out <- .Fortran("mweave", as.double(x), as.double(dw),
                    as.double(swz), as.integer(ndata),
                    as.integer(nord), as.double(ssqswz),
                    cntr=double(1), as.double(dt_),
                    spz=double(1), varc=double(1),
                    PACKAGE='multitaper')
    return(list(cntr=out$cntr, spz=out$spz, varc=out$varc))
}

.HF4mp1 <- function(cft, swz, nord, ssqswz) {
    
    ## vectorized when there are none to skip... 
    ## Dave's [82] mu hat of f in eqn (13.5)
    ## cmv   Complex Mean Value,
    cmv <- (cft %*% swz) /ssqswz
    ssqave <-  abs(cmv)**2*ssqswz
    ##cftEst <- t(swz %*% t(cmv))
    swz <- as.matrix(swz)
    ##removed extra crossproduct 
    ##cftEst <- t(tcrossprod(swz,cmv))
    
    ssqres <- apply( abs(cft - (cmv %*% t(swz)))**2,
                    1, sum)
    F_<- Re((nord-1)*ssqave/ssqres)
    
    return(list(Ftest=F_,cmv=cmv))
}


.mw2wta <- function(sa, nfreq, nord,
                    var, dt_, ev, evp=(1-ev),
                    tol=.03, maxadaptiveiteration=100) {

    out <- .Fortran("mw2wta", as.double(sa),
                    wt=matrix(as.double(0), nfreq, nord),
                    as.integer(nfreq), as.integer(nord),
                    s=double(nfreq), as.double(ev), as.double(evp),
                    dofs=double(nfreq), dofav=double(1),
                    as.double(var), as.double(dt_),
                    as.double(tol),
                    as.integer(maxadaptiveiteration),
                    mxiter=integer(1), aviter=double(1),
                    PACKAGE='multitaper')
    
    return(list(s=out$s, wt=out$wt, dofs=out$dofs, dofav=out$dofav,
                mxiter=out$mxiter, aviter=out$aviter))
}

.mw2jkw <-  function(sa, nfreq, nord, var, dt_, ev,
                     evp=(1-ev), tol=.03,
                     maxadaptiveiteration=100) {
    
    nordP2 <-  nord+2
    out <- .Fortran("mw2jkw", as.double(sa),
                    wt=matrix(as.double(0), nfreq, nord),
                    as.integer(nfreq), as.integer(nord),
                    s=double(nfreq), as.double(ev),
                    as.double(evp), dofs=double(nfreq),
                    dofav=double(1), as.double(var),
                    as.double(dt_), as.double(tol),
                    sjk=double(nordP2), varjk=double(nfreq),
                    bcjk=double(nfreq),
                    matrix(as.double(0), nord, nordP2),
                    double(nordP2), double(nord),
                    as.integer(maxadaptiveiteration),
                    mxiter=integer(1),
                    PACKAGE='multitaper')
                   
    
    return(list(s=out$s, wt=out$wt, dofs=out$dofs, dofav=out$dofav,
                mxiter=out$mxiter,
                varjk=out$varjk, bcjk=out$bcjk, sjk=out$sjk))
}
