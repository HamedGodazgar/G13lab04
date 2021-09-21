linreg <-
function(formula , data) {
  
  # use model.matrix() to create the matrix X (independent variables)
  X <- model.matrix(formula,data)
  
  # pick out the dependent variable y using all.vars()
  y <- data[all.vars(formula)[1]]
  
  ## QR factorization
  QR <- qr(X)
  r <- QR$rank
  piv <- QR$pivot[1:r]
  
  ## estimate identifiable coefficients
  betahat <- solve.qr(QR, y)
  betahat
  ## fitted values
  yhat <- as.vector(c(X[, piv] %*% betahat))
  y <-as.matrix( data[all.vars(formula)[1]])
  ## residuals
  resi <- as.vector(y - yhat)
  
  ## degree of freedom
  dof <- nrow(X) - r
  
  ## error variance
  se2 <- c(crossprod(resi)) / dof
  
  ## variance-covariance for coefficients
  v <- chol2inv(QR$qr, r) * se2
  
  ## t-value
  tval <- betahat / sqrt(diag(v))
  yy <- data[all.vars(formula)[1]]
  
  ## return
  return_obj <- list(Coefficients = betahat, fitted.values = yhat,
                     residuals = resi, degfreedom = dof,
                     vcov = v, varr = sqrt(se2), t_value = tval,call=match.call(),y=yy)
  
  class(return_obj) <- 'linreg'
  return(return_obj)
  
}
