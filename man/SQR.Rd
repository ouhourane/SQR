\name{SQR}
\alias{SQR}
\title{
Fits the regularization paths for the penalized quantile regression
}
\description{
Fits a regularization path for the penalized quantile regression at a sequence of regularization parameters lambda.
}
\usage{
SQR(x, y, nlambda = 100, check = c("f1", "f2"), method = c("Lasso", "Scad", "Mcp"), lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04), lambda = NULL, lambda2 = 0, pf = rep(1, nvars), pf2 = rep(1, nvars), exclude, dfmax = nvars + 1, pmax = min(dfmax * 1.2, nvars), standardize = TRUE, eps = 1e-06, maxit = 1e+06, kk = 1, c = 1, taux = 0.5, gamma = 3)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{x}{matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}
\item{y}{response variable: vector of N observations.}
\item{nlambda}{the number of \code{lambda} values - default is 100.}

\item{check}{a character string specifying the ckeck function approximation to use, valid options are:
		\itemize{
		\item \code{"f1"}: the check function 1,
		\item \code{"f2"}: the check function 2.
		} Default is \code{"f1"}.}

\item{method}{a character string specifying the function penality to use, valid options are:
		\itemize{
		\item \code{"Lasso"}: the Lasso penality,
		\item \code{"Scad"}: the SCAD penality,
                \item \code{"Mcp"}: the MCP penality,
		} Default is \code{"Lasso"}.}

\item{lambda.factor}{The factor for getting the minimal lambda in \code{lambda} sequence, where \code{min(lambda)} = \code{lambda.factor} * \code{max(lambda)}.  \code{max(lambda)} is the smallest value of \code{lambda} for which all coefficients are zero. The default depends on the relationship between \eqn{N} (the number of rows in the matrix of predictors) and \eqn{p}{p} (the number of predictors). If \eqn{N > p}, the default is \code{0.0001},
		close to zero.  If \eqn{N<p}, the default is \code{0.01}.
		A very small value of \code{lambda.factor} will lead to a saturated fit. It takes no effect if there is user-defined \code{lambda} sequence.} 
		\item{lambda}{a user supplied \code{lambda} sequence. Typically, by leaving this option unspecified users can have 
		the program compute its own \code{lambda} sequence based on
		\code{nlambda} and \code{lambda.factor}. Supplying a value of
		\code{lambda} overrides this. It is better to supply
		a decreasing sequence of \code{lambda} values than a single (small) value, if not, the program will sort user-defined \code{lambda} sequence in decreasing order automatically.}

\item{lambda}{
regularization parameter \eqn{\lambda_1}{lambda1} for the l1 penalty of the 
coefficients
}

\item{lambda2}{regularization parameter \eqn{\lambda_2}{lambda2} for the quadratic penalty of the 
coefficients.}

\item{pf}{L1 penalty factor of length \eqn{p}{p} used for adaptive LASSO or adaptive elastic net. Separate L1 penalty weights can be applied to each coefficient of \eqn{\beta}{beta} to allow
		differential L1 shrinkage. Can be 0 for some variables, which implies
		no L1 shrinkage, and results in that variable always being included in the
		model. Default is 1 for all variables (and implicitly infinity for
		variables listed in \code{exclude}).}
		
\item{pf2}{L2 penalty factor of length \eqn{p}{p} used for adaptive LASSO or adaptive elastic net. Separate L2 penalty weights can be applied to each coefficient of \eqn{\beta}{beta} to allow
		differential L2 shrinkage. Can be 0 for some variables, which implies
		no L2 shrinkage. Default is 1 for all variables.}

\item{exclude}{indices of variables to be excluded from the
		model. Default is none. Equivalent to an infinite penalty factor.}

\item{dfmax}{limit the maximum number of variables in the
		model. Useful for very large \eqn{p}, if a partial path is desired. Default is \eqn{p+1}.}

\item{eps}{convergence threshold for coordinate majorization descent. Each inner
		coordinate majorization descent loop continues until the relative change in any
		coefficient (i.e. \eqn{\max_j(\beta_j^{new}-\beta_j^{old})^2}{max(j)(beta_new[j]-beta_old[j])^2}) is less than \code{eps}. For HHSVM i.e. \code{method="hhsvm"}, it is \eqn{\frac{2}{\delta}\max_j(\beta_j^{new}-\beta_j^{old})^2}{2*max(j)(beta_new[j]-beta_old[j])^2/delta}. Defaults value is \code{1e-8}.}

\item{maxit}{maximum number of outer-loop iterations allowed at fixed lambda value. Default is 1e6. If models do not converge, consider increasing \code{maxit}.}

\item{kk}{
the pamameter k in the appriximation check function 2
}
  \item{c}{
the pamameter c in the appriximation check function 1
}

  \item{taux}{
   quantile to be targeted. Must be in (0,1).
}

\item{gamma}{
regularization parameter for the Scad or Mcp penalty of the 
coefficients
}

}
\details{
Note that the objective function in \code{SQR} is
\deqn{Loss(y, X, beta))/N + lambda1 * |beta| + 0.5 * lambda2 * beta^2}
where the penalty is a combination of L1 and L2 term.
}

\value{
An object with S3 class \code{\link{SQR}}.
		\item{call}{the call that produced this object}
		\item{b0}{intercept sequence of length \code{length(lambda)}}
		\item{beta}{a \code{p*length(lambda)} matrix of coefficients, stored as a sparse matrix (\code{dgCMatrix} class, the standard class for sparse numeric matrices in the \code{Matrix} package.). To convert it into normal type matrix use \code{as.matrix()}.}
		\item{lambda}{the actual sequence of \code{lambda} values used}
		\item{df}{the number of nonzero coefficients for each value of
		\code{lambda}.}
		\item{dim}{dimension of coefficient matrix (ices)}
		\item{npasses}{total number of iterations (the most inner loop) summed over all lambda values}
		\item{jerr}{error flag, for warnings and errors, 0 if no error.}
}
\author{Abdallah Mkhadri and Mohamed Ouhourane
Maintainer: Abdallah Mkhadri <mkhadri@uca.ma> Mohamed OUHOURANE <mohamed.ouhourane@gmail.com> Karim Oulkacha <oualkacha.karim@uqam.ca> }
\references{
Abdallah Mkhadri, Mohamed Ouhourane, Karim Oulkacha: "A coordinate descent algorithm for computing the penalized smooth quantile regression" 
}
\examples{
x <- matrix(rnorm(100*10),100,10)
y <- 2*x[,1]-x[,2]+rnorm(100,0,1)
fit <- SQR(x=x,y=y, check = "f2", method = "Lasso")
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

