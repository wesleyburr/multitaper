%     The multitaper R package
%     Multitaper and spectral analysis package for R
%     Copyright (C) 2011 Karim Rahim 

%     This file is part of the multitaper package for R.

%     The multitaper package is free software: you can redistribute it and
%     or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 2 of the License, or
%     any later version.

%     The multitaper package is distributed in the hope that it will be 
%     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
%     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.

%     You should have received a copy of the GNU General Public License
%     along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

%     If you wish to report bugs please contact the author. 
%     karim.rahim@gmail.com
%     112 Jeffery Hall, Queen's University, Kingston Ontario
%     Canada, K7L 3N6

\name{spec.mtm}
\alias{spec.mtm}
\title{Compute and plot the multitaper spectrum estimate}
\description{
   Computes and plots the adaptive multitaper spectrum estimate.
}
\usage{
spec.mtm(timeSeries, k, nw, nFFT = "default",
centreWithSlepians = TRUE, dpssIN = NULL, returnZeroFreq = TRUE, Ftest =
FALSE, jackknife = FALSE, jkCIProb = 0.95, maxAdaptiveIterations = 100,
plot = TRUE, na.action = na.fail, returnInternals = FALSE, ...) 
}
\arguments{
  \item{timeSeries}{A time series of equally spaced data, this can be
created by the ts() function where deltat is specified.}
  \item{k}{k a positive integer, the number of tapers, often 2*nw.}
  \item{nw}{nw a positive double precision number, the time-bandwidth
    parameter.}
  \item{nFFT}{This function pads the data before computing the fft. nFFT
indicates the total length of the data after padding.}
  \item{centreWithSlepians}{The data series is centred using an
    expansion into the Slepian's (dpss).}
  \item{dpssIN}{Allows the user to enter a dpss object which has already
    been created. This can save computation time when dpss sequences wth
 the same bandwidth parameter and same number of tapers are used repeatedly.}
  \item{returnZeroFreq}{Boolean variable indicating if the zeroth
frequency (DC component) should be returned.}
  \item{Ftest}{Boolean variable indicating if the Ftest result should be
returned.}
  \item{jackknife}{Boolean variable indicating if jackknifed confidence
intervals should be returned.}
  \item{jkCIProb}{Decimal value indicating the jackknife probability to
use in determining jackknife confidence intervals, the default returns a
95\% confidence interval.}
  \item{maxAdaptiveIterations}{Maximum number of iterations in the
adaptive multitaper calculation, generally convergence is quick.}
  \item{plot}{Boolean variable indicating if the spectrum should be plot}
  \item{na.action}{Action to take if NA's exist in the data, the default
    is to fail.}
  \item{returnInternals}{Aause the function to return the
    eigencoefficients, and/or weighted eigencoefficients. For
    example mtm.coh requires mtm objects with internals.}
  \item{...}{Additional parameters, such as xaxs="i" which are passed to
    the plotting function.} 
}
\references{
Thomson, D.J (1982)
Spectrum estimation and harmonic analysis. \emph{Proceedings of the IEEE}
Volume \bold{70}, number 9, pp. 1055--1096.

Percival, D.B. and Walden, A.T. (1993)
\emph{Spectral analysis for physical applications}
Cambridge University Press. 
}
\examples{
data(willamette)
spec.mtm(willamette, nw=4, k=8)
spec.mtm(willamette, nw=4, k=8, nFFT=2048)

}
\keyword{multitaper}
