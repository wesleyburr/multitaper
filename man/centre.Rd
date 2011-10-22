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

\name{centre}
\alias{centre}
\title{Centres the data or time series.}
\description{
  Centres the data using the Slepian sequences if the bandwidth
  parameters, nw, and number of tapers, k, is specified, otherwise
  subtracts the mean or trim mean.
}
\usage{
centre(x, nw = NULL, k = NULL, deltaT = NULL, trim = 0)

}
\arguments{
  \item{x}{x the data as a vector or as a time series.}
  \item{nw}{nw the Slepian bandwidth parameter.}
  \item{k}{k the number of tapers used, often 2*nw.}
  \item{deltaT}{deltaT required if the data is a vector and not a time
    series if the Slepians are used to centre the data.}
  \item{trim}{trim [only used if nw and k are not specified]
    the fraction (0 to 0.5) of observations to be trimmed from
    each end of `x' before the mean is computed.  Values of trim
    outside that range are taken as the nearest endpoint.}
}
\references{
  Thomson, D.J (1982)
  Spectrum estimation and harmonic analysis. \emph{Proceedings of the IEEE}
  Volume \bold{70}, number 9, pp. 1055--1096.
  
  Slepian, D. (1978)
  Prolate spheroidal wave functions, Fourier analysis, and
  uncertainty. V--The discrete case. \emph{Bell System Technical Journal}
  Volume \bold{57}, pp. 1371--1430.
}
\examples{
data(willamette)
centredWDat <- centre(willamette, nw=4, k=8, deltaT=1)
}
\keyword{math}
