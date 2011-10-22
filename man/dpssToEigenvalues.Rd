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

\name{dpssToEigenvalues}
\alias{dpssToEigenvalues}
\title{Compute eigenvalues for the Discrete Prolate Spheroidal Sequences
(dpss)}
\description{
   Compute eigenvalues for the Discrete Prolate Spheroidal Sequences. The
   method used here is described in chapter 8 of Percival and Walden (1993).
}
\usage{
dpssToEigenvalues(v, nw)
}
\arguments{
  \item{v}{v a matrix of dpss, with each column representing a dpss of a
    different order.}
  \item{nw}{nw a positive double precision number, the time-bandwidth
    parameter.}
}
\references{
Percival, D.B. and Walden, A.T. (1993)
\emph{Spectral analysis for physical applications.}
Cambridge University Press.
}
\examples{
## dpss(10,4,4)
dpss1 <- dpss(10,4,4, returnEigenvalues=FALSE)$v
dpssToEigenvalues(dpss1,4)
}
\keyword{math}
