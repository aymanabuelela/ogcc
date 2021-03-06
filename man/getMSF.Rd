% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ogcc.R
\name{getMSF}
\alias{getMSF}
\title{Get minimum set of features (getMSF)}
\usage{
getMSF(model = "types", d = "RSEM")
}
\arguments{
\item{model}{is a prediction model name from \code{\link{ogcc}}. See \code{model} argument.}

\item{d}{a string indicating the RNA-Seq measurement type; either 'RSEM' or 'RPKM'. 'RSEM' by default.}
}
\value{
a character vector of the required features for a working model.
}
\description{
\code{getMSF} prints the minimum set of features required by a model from \code{ogcc}.
}
\examples{
msf <- getMSF("normal_tumor")
print(msf)

}
\author{
Ayman Abuelela; ayman.elkhodiery@kaust.edu.sa
}
\seealso{
\code{\link{ogcc}} and \code{\link{getLabels}}
}

