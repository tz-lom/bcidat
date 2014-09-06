\name{bcidat-package}
\alias{bcidat-package}
\alias{bcidat}
\docType{package}
\title{
Adds support of BCI2000 .dat files
}
\description{
Based on top of BCI2000 source, this packages allows reading of BCI2000 .dat files that contains unprocessed EEG recordings.
}
\details{
\tabular{ll}{
Package: \tab bcidat\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2014-09-04\cr
License: \tab GPL-3\cr
}

Package allows reading of EEG data from BCI2000 .dat files.
}
\author{
Nuzhdin Yury

Maintainer: Nuzhdin Yury <nuzhdin.urii@gmail.com>
}
\references{
Visit \href{http://www.bci2000.org/wiki/index.php/Main_Page}{BCI2000 Wiki} for more information about BCI2000
}
\keyword{ package }
\examples{
data <- load_bcidat('filename')
}
