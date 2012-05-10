%     The multitaper R package
%     Multitaper and spectral analysis package for R
%     Copyright (C) 2011 Karim Rahim 
%
%     Written by Karim Rahim and Wesley Burr.
%
%     This file is part of the multitaper package for R.
%
%     The multitaper package is free software: you can redistribute it and
%     or modify it under the terms of the GNU General Public License as 
%     published by the Free Software Foundation, either version 2 of the 
%     License, or any later version.
%
%     The multitaper package is distributed in the hope that it will be 
%     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
%     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with multitaper.  If not, see <http://www.gnu.org/licenses/>.
%
%     If you wish to report bugs please contact the author. 
%     karim.rahim@gmail.com
%     112 Jeffery Hall, Queen's University, Kingston Ontario
%     Canada, K7L 3N6

\name{prewh}
\title{Prewhiten a given time series}
\description{
    Computes a robust estimate of the autocorrelation of a time series
    using the inverse multitaper spectrum estimate technique. Using a 
    user-parameter for maxlag, prewhitens a given time series.
}
\usage{
prewh(x,dT,maxlag)
}
\arguments{
  \item{x}{Time series, required to be contiguous. Typically a numeric vector,
          not a ts object.}
  \item{dT}{Time step for series.}
  \item{maxlag}{Integer giving maximum number of lags to be used in the
        prewhitening.}
}
\examples{
xd <- rnorm(1000)
yd <- prewh(x=xd,dT=1,maxlag=10)
}
