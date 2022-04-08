\name{coef.SQR}
\alias{coef.SQR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{get coefficients or make coefficient predictions from an "SQR" object.}

\description{
Computes the coefficients or returns a list of the indices of the nonzero coefficients at the requested values for \code{lambda} from a fitted \code{\link{SQR}} object.
}

\arguments{
  	        \item{object}{fitted \code{\link{SQR}} model object.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the entire sequence used to
		create the model.}
		\item{type}{type \code{"coefficients"} computes the coefficients at the requested
		values for \code{s}. Type}
		\item{\dots}{not used. Other arguments to predict. } 
}
\details{\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{coef} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of coefficients from both left and right \code{lambda} indices.}

\author{Abdallah Mkhadri and Mohamed Ouhourane
Maintainer: Abdallah Mkhadri <mkhadri@uca.ma> Mohamed OUHOURANE <mohamed.ouhourane@gmail.com> Karim Oulkacha <oualkacha.karim@uqam.ca> }
\references{
Abdallah Mkhadri, Mohamed Ouhourane, "A coordinate descent algorithm for computing the penalized smooth quantile regression" 
}
\examples{
x <- matrix(rnorm(100*10),100,10)
y <- 2*x[,1]-x[,2]+rnorm(100,0,1)
fit<- SQR(x=x,y=y,check="f1" method = "Lasso")
coef(fit,type="coef",s=c(0.1,0.005))
coef(fit,type="nonzero")
}

\keyword{models}
\keyword{regression}
\keyword{elastic net}
\keyword{lasso}
\keyword{Scad}
\keyword{MCP}
\keyword{penalized quantile regression}
\keyword{coordinate descente algorithm}
