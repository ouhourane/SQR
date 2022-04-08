\name{plot.SQR}
\alias{plot.SQR}
\title{Plot coefficients from a "SQR" object}

\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted \code{\link{SQR}} object. This function is modified based on the \code{plot} function from the \code{gcdnet} package.
}
\usage{
plot.SQR(x, xvar = c("norm", "lambda"), color = FALSE, label = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  	\item{x}{fitted \code{\link{SQR}} model}
		\item{xvar}{what is on the X-axis. \code{"norm"} plots against the
		L1-norm of the coefficients, \code{"lambda"} against the log-lambda
		sequence.}
		\item{color}{if \code{TRUE}, plot the curves with rainbow colors. \code{FALSE} is gray colors. Default is \code{FALSE}}
		\item{label}{if \code{TRUE}, label the curves with variable sequence numbers. Default is \code{FALSE}}
		\item{\dots}{other graphical parameters to plot}
}

\details{
A coefficient profile plot is produced.
}
\author{Abdallah Mkhadri and Mohamed Ouhourane
Maintainer: Abdallah Mkhadri <mkhadri@uca.ma> Mohamed OUHOURANE <mohamed.ouhourane@gmail.com> Karim Oulkacha <oualkacha.karim@uqam.ca> }
\references{
Abdallah Mkhadri, Mohamed Ouhourane, Karim Oulkacha "A coordinate descent algorithm for computing the penalized smooth quantile regression" 
}
\examples{
x <- matrix(rnorm(100*10),100,10)
y <- 2*x[,1]-x[,2]+rnorm(100,0,1)
fit <- SQR(x=x,y=y, check = "f1", method = "Lasso")
plot(fit) # plots against the L1-norm of the coefficients
}

\keyword{models}
\keyword{regression}
\keyword{elastic net}
\keyword{lasso}
\keyword{Scad}
\keyword{MCP}
\keyword{penalized quantile regression}
\keyword{coordinate descente algorithm}
