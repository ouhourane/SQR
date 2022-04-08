qgcdpathMcp1 <- function(x, y, nlam, flmin, ulam, isd, 
                          eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, kk, taux, nobs, nvars, 
                          vnames, gamma) {
  #################################################################################
  #data setup
  y <- as.numeric(y)
  if (kk < 0) 
    stop("c must be non-negative")
  kk <- as.double(kk)
  if (taux > 1) 
    stop("taux must be betwwen 0 and 1")
  if (taux < 0) 
    stop("taux must be betwwen 0 and 1")
  taux <- as.double(taux)
  gamma <- as.double(gamma)
  #################################################################################
  # call Fortran core
  fit <- .Fortran("quantileMcpNET1", gamma, kk, taux, lam2, nobs, nvars, 
                  as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                  flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                  beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                  alam = double(nlam), npass = integer(1), jerr = integer(1))
  #################################################################################
  # output
  outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
  outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
  class(outlist) <- c("cdqpath")
  outlist
} 
