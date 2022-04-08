cv.SQR <-
function (x, y, K = 5,lambda2=0, kk=1, c=1, taux=0.5, check = c("f1","f2"),
                     method = c("Lasso","Scad","Mcp"),plot.it=T)
  {
  
    check <- match.arg(check)
    method <- match.arg(method)
    n=dim(x)[1]
    p=dim(x)[2]
    all.folds <- split(sample(1:length(y)), rep(1:K, length = length(y)))
    object=SQR(x,y=y,lambda2=lambda2,check=check,c=c,taux=taux,method = method)
    lambda <- object$lambda
    residmat <- matrix(0, length(lambda), K)
    for (i in seq(K)) {
      omit <- all.folds[[i]]
      fit <- SQR(x=x[-omit, ],y=y[-omit],lambda2=lambda2,lambda=lambda,check=check,method = method,kk=kk, c=c,taux=taux)
      beta0 <- as.matrix(fit$b0)
      beta <- fit$beta
      ypred <- t((x[omit, ]%*%beta))+beta0
      residmat[, i] <- apply((y[omit] - t(ypred))^2, 2, mean)
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    finalfit <- object

    if (plot.it) {
      matplot(lambda, cv, type = "l",main = paste("CV",method) ,ylim = range(cv, cv + cv.error, cv - cv.error))
        error.bars(lambda, cv + cv.error, cv - cv.error, width = 1/length(lambda),xlab = "lambda")
    }
    
    list(cv=cv,cv.error=cv.error,finalfit=finalfit,all.folds=all.folds,lambda = object$lambda)
    }
