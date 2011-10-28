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
%     along with multitaper.  If not, see <http://www.gnu.org/licenses/>.

%     If you wish to report bugs please contact the author. 
%     karim.rahim@gmail.com
%     112 Jeffery Hall, Queen's University, Kingston Ontario
%     Canada, K7L 3N6

\name{dropFreqs}
\alias{dropFreqs}
\alias{dropFreqs.default}
\alias{dropFreqs.spec}
\alias{dropFreqs.mtm}
\alias{dropFreqs.mtm.coh}
\title{Compute and plot the multitaper spectrum estimate}
\description{
   A utility function to reduce the frequencies in a spectral
   estimate. This utility is used before calling plot, to increase the
   frequency resolution of a plot by dropping frequencies outside a
   particular band of interest. This is not a filter, it is simply a
   uitilty to allow R to restrict a plot to a certain band.
}
\usage{
dropFreqs(spec, minFreq, maxFreq)
}
\arguments{
  \item{spec}{spec a spectrum estimate returned, of class spec, mtm, or mtm.coh.}
  \item{minFreq}{minFreq the minium frequency to keep, all frequency,
    and spectral values below this will be dropped.}
  \item{maxFreq}{maxFreq the maximum frequency to keep, all frequency,
    and spectral values below this will be dropped.}
}
\examples{
data(willamette)
mtm1 <- spec.mtm(willamette, nw=4, k=8, plot=FALSE)
mtm2 <- dropFreqs(mtm1, .1, .4)
plot(mtm2)
}
\keyword{multitaper}
