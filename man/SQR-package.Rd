\name{SQR-package}
\alias{SQR-package}
\alias{SQR}
\docType{package}
\title{
A coordinate descent algorithm for computing the penalized
smooth quantile regression
}
\description{
This package implement a coordinate descent algorithm for
computing the penalized smooth quantile regression (cdaSQR). The main idea is to approx-
imate the objective check function, which is not differentiable at zero, with a modified check
function which is differentiable at zero. Then, using the maximization-minimization trick of
the gcdnet algorithm, we update each coefficient simply and efficiently.
}
\author{Abdallah Mkhadri and Mohamed Ouhourane
Maintainer: Abdallah Mkhadri <mkhadri@uca.ma> Mohamed OUHOURANE <mohamed.ouhourane@gmail.com> Karim Oulkacha <oualkacha.karim@uqam.ca> }
\references{
Abdallah Mkhadri, Mohamed Ouhourane, Karim Oulkacha: "A coordinate descent algorithm for computing the penalized smooth quantile regression" 
}
\examples{
x <- matrix(rnorm(100*10),100,10)
y <- 2*x[,1]-x[,2]+rnorm(100,0,1)
fit <- SQR(x=x,y=y, check = "f2", method = "Scad")
plot(fit)
}

\keyword{models}
\keyword{regression}
\keyword{elastic net}
\keyword{lasso}
\keyword{Scad}
\keyword{MCP}
\keyword{penalized quantile regression}
\keyword{coordinate descente algorithm}
