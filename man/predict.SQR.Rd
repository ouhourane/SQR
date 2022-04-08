\name{predict.SQR}

\alias{predict.SQR}

\title{make predictions from a "SQR" object.}

\description{
Similar to other predict methods, this functions predicts fitted values and class labels from a fitted \code{SQR} object.
}

\usage{
predict.SQR(object, newx, s = NULL)
}

\arguments{
    \item{object}{fitted \code{\link{SQR}} model object.}
		\item{newx}{matrix of new values for \code{x} at which predictions are
		to be made. NOTE: \code{newx} must be a matrix, \code{predict} function does not accept a vector or other formats of code{newx}.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the entire sequence used to
		create the model. }
}

\details{\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{predict} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of predicted values from both left and right \code{lambda} indices. }
\author{Abdallah Mkhadri and Mohamed Ouhourane
Maintainer: Abdallah Mkhadri <mkhadri@uca.ma> Mohamed OUHOURANE <mohamed.ouhourane@gmail.com> Karim Oulkacha <oualkacha.karim@uqam.ca> }
\references{
Abdallah Mkhadri, Mohamed Ouhourane, Karim Oulkacha: "A coordinate descent algorithm for computing the penalized smooth quantile regression" 
}
\examples{
x <- matrix(rnorm(100*10),100,10)
y <- 2*x[,1]-x[,2]+rnorm(100,0,1)
fit=SQR(x=x,y=y, check = "f2", method = "Lasso")
predict(fit,newx=x[2:10,])
}

\keyword{models}
\keyword{regression}
\keyword{elastic net}
\keyword{lasso}
\keyword{Scad}
\keyword{MCP}
\keyword{penalized quantile regression}
\keyword{coordinate descente algorithm}
