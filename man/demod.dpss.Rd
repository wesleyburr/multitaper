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

\name{demod.dpss}
\alias{demod.dpss}
\title{Computes complex demodulates using multiple taper techniques}
\description{
  Computes complex demodulate of a given series around a given central
  frequency using multiple taper techniques. Returns amplitude, phase,
  and complex demodulate.
}
\usage{
demod.dpss(x,centreFreq,NW,blockLen,stepSize=1,wrapphase=TRUE) 
}
\arguments{
  \item{x}{Time series, required to be contiguous.}
  \item{centreFreq}{Frequency around which to demodulate.}
  \item{NW}{Parameter controlling time-bandwidth.}
  \item{blockLen}{Length of sub-block to use; demodulate is computed on 
  each block in turn.}
  \item{stepSize}{If true the appropriate eigenvalues are
    calculated and returned using the function dpssToEigenvalues.}
  \item{wrapphase}{If true, routine wraps phases around +/-360 degree 
  boundaries.}
}
\references{
}
\examples{
}
\keyword{math}
