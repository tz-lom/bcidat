\name{load_bcidat}
\alias{load_bcidat}
\docType{package}
\title{
Reads .dat file
}
\description{
Loads signal, state and parameters from .dat file
}
\usage{
load_bcidat(file, raw = FALSE)	
}
\arguments{
  \item{file}{
    Name of file to be loaded.
    Extension can be omitted, so function will try to open `file.dat` if `file` is not existing. 
  }
  \item{raw}{
    Whether load raw data, or calibrated. 
  }
}
\value{
  \item{signal}{
    Matrix of the dimension samples*channels with EEG data
  }
  \item{states}{
    Matrix with state values. Number of rows corresponds to number of samples in signal.
  }
  \item{parameters}{
    List of parameters. Values can be characters, matrices of characters or lists of lists of anything else.
  }
}
\examples{
\dontrun{
data <- load_bcidat('record.dat')
}
}
