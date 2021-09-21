

#' @author Hamed , Omid
#' @description Implementing Liear Regression using QR Decomposition
#' @examples 
#' mod_object <- linreg(Petal.Length~Species, data = iris)
#' print(mod_object)
#' resid(mod_object)
#' @export linreg
#' @exportS3Method 
#' print()
#' plot()
#' resid()
#' pred()
#' coef()
#' summary()
#' @import ggplot2
#' @name linreg
#' @param formula
#' @param data
#' @references \url{http://staff.www.ltu.se/~jove/courses/c0002m/least_squares.pdf}
#' @title linreg
#' @usage linreg(formula,data)

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




print.linreg <- function(obj){
  cat("Call:\n")
  print(obj$call)
  cat("\nCoefficients:\n")
  print(obj$Coefficients)
  
}
print(mod_object)

resid.linreg <- function(obj){
  return(obj$residuals)
  
}
resid(mod_object)
pred <- function(obj){
  return(obj$fitted.values)
  
}
pred(mod_object)

coef.linreg <- function(obj){
  names(obj$Coefficients)
  
  
}
coef(mod_object)




summary.linreg <- function(obj)
{
  se <- sqrt(diag(obj$vcov))
  tval <- obj$t_value
  TAB <- cbind(Estimate = obj$Coefficients,
               StdErr = se,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=obj$degfreedom))
  res <- list(call=obj$call,
              coefficients=TAB)
  class(res) <- "summary.linreg"
  print(res)
  cat("\nResidual standard error:" , obj$varr ,"on", obj$degfreedom ,"degrees of freedom")
}
summary(mod_object)







plot.linreg <- function(obj) {
  library(ggplot2)
  par(mfrow=c(2,1))
  ggplot(obj$Petal.Length,aes(x = obj$fitted.values, y = obj$residuals),
         color="steelblue")+
    geom_point(shape=1,size=4) + geom_smooth(method="lm", colour="red",
                                             se = FALSE)+
    labs(title="Residuals vs Fitted",x="Fitted values
      linreg(Petal.Length~Species, data = iris)", y="Residuals")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "white", colour = "grey50"))
  
  
  
  
}
plot(mod_object)