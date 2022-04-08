\name{cv.SQR}
\alias{cv.SQR}
\title{Cross-validation for SQR}
\description{Does K-fold cross-validation for SQR, produces a plot,
and returns a value for \code{lambda}. This function is modified based on the \code{cv} function from the \code{VISA}.}

\usage{
cv.SQR(x, y, K = 5, lambda2 = 0, kk = 1, c = 1, taux = 0.5, check = c("f1", "f2"), plot.it = T)
}
\arguments{
\item{x}{\code{x} matrix as in \code{\link{SQR}}.}
\item{y}{response}
\item{K}{number of folds - default is 5. Although \code{K=5}
		can be as large as the sample size (leave-one-out CV), it is not
		recommended for large datasets. Smallest value allowable is \code{nfolds=3}.}
\item{lambda2}{regularization parameter \eqn{\lambda_2}{lambda2} for the quadratic penalty of the  coefficients.}
\item{kk}{
the pamameter k in the approximation check function 2
}
\item{c}{
the pamameter c in the approximation check function 1
}
\item{taux}{
   quantile to be targeted. Must be in (0,1).
}
\item{check}{a character string specifying the ckeck function approximation to use, valid options are:
		\itemize{
		\item {\code{"f1"}: The check function 1}
		\item {\code{"f2"}: The check function 2}
		} 
}
\item{plot.it}{ Plot it?}
}


\details{The function runs \code{\link{SQR}}. the first to get the \code{lambda} sequence, and then the remainder to compute the fit with each of the folds omitted. The average error and standard deviation over the folds are computed. 
}

\value{an object of class \code{cv.SQR} is returned, which is a
list with the ingredients of the cross-validation fit.

		\item{cv}{the mean cross-validated error - a vector of length}
		\item{cv.error}{estimate of standard error of \code{cv}.}
		\item{finalfit}{An object with S3 class \code{SQR.}}
		\item{all.folds}{ The folds used in cross-validation.}
		\item{lambda}{the values of \code{lambda} used in the fits.}
}
\author{Abdallah Mkhadri and Mohamed Ouhourane
Maintainer: Abdallah Mkhadri <mkhadri@uca.ma> Mohamed OUHOURANE <mohamed.ouhourane@gmail.com> Karim Oulkacha <oualkacha.karim@uqam.ca> }
\references{
Abdallah Mkhadri, Mohamed Ouhourane, Karim Oulkacha: "A coordinate descent algorithm for computing the penalized smooth quantile regression" 
}
\examples{
x <- matrix(rnorm(100*10),100,10)
y <- 2*x[,1]-x[,2]+rnorm(100,0,1)
cv=cv.SQR(x=x,y=y, check = "f1", method = "Lasso")
}

\keyword{models}
\keyword{regression}
\keyword{elastic net}
\keyword{lasso}
\keyword{Scad}
\keyword{MCP}
\keyword{penalized quantile regression}
\keyword{coordinate descente algorithm}

