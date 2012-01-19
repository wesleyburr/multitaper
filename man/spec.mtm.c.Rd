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

\name{spec.mtm.c}
\alias{spec.mtm.c}
\title{Compute and plot (complex data) multitaper spectrum estimates}
\description{
   Computes and plots adaptive or nonadaptive multitaper spectrum estimates from
   contiguous *complex* time series objects
}
\usage{
spec.mtm.c(timeSeries, nw=4.0, k=8, nFFT="default", plot=TRUE, na.action=na.fail,
                     returnInternals=TRUE, dtUnits=c("default"), dT=NULL, ...) 
}
\arguments{
  \item{timeSeries}{A complex time series of equally spaced data, this can be
created by the ts() function where deltat is specified.}
  \item{nw}{nw a positive double precision number, the time-bandwidth
    parameter.}
  \item{k}{k a positive integer, the number of tapers, often 2*nw.}
  \item{nFFT}{This function pads the data before computing the fft. nFFT
indicates the total length of the data after padding.}
  \item{plot}{Boolean variable indicating if the spectrum should be plotted.}
  \item{na.action}{Action to take if NAs exist in the data, the default
    is to fail.}
  \item{returnInternals}{Return the weighted eigencoefficients, complex mean values, 
    and so on. These are necessary for extensions to the multitaper, including 
    magnitude-squared coherence (function mtm.coh in this package).}
  \item{dtUnits}{Allows indication of the units of delta-t for accurate frequency axis
  labels.}
  \item{dT}{Time step for observations. If not in seconds, dtUnits should be set to indicate
  the proper units for plot labels.}
  \item{...}{Additional parameters, such as xaxs="i" which are passed to
    the plotting function. Not all parameters are supported.} 
}
\references{
Thomson, D.J (1982)
Spectrum estimation and harmonic analysis. \emph{Proceedings of the IEEE}
Volume \bold{70}, Number 9, pp. 1055--1096.

Percival, D.B. and Walden, A.T. (1993)
\emph{Spectral analysis for physical applications}
Cambridge University Press. 
}
\examples{
}
\keyword{multitaper}
