

#' @author Hamed , Omid
#' @description Implementing Liear Regression using QR Decomposition
#' @examples 
#' mod_object <- linreg(Petal.Length~Species, data = iris)
#' print(mod_object)
#' resid(mod_object)
#' @export linreg
#' @import ggplot2
#' @import cowplot
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
  std <- sqrt(abs(resi))
  ## return
  return_obj <- list(Coefficients = betahat, fitted.values = yhat,
                     residuals = resi, degfreedom = dof,
                     vcov = v, varr = sqrt(se2), t_value = tval,
                     call=match.call(), y=yy, standardized= std)
  
  class(return_obj) <- 'linreg'
  return(return_obj)
  
}

mod_object <- linreg(Petal.Length~Species, data=iris)
linreg_mod <- linreg(Petal.Length~Sepal.Width+Sepal.Length, data=iris)

# print method
print.linreg <- function(obj){
  cat("Call:\n")
  print(obj$call)
  cat("\nCoefficients:\n")
  print(obj$Coefficients)
  
}

# resid method
resid.linreg <- function(obj){
  return(obj$residuals)
  
}

# pred method
pred <- function(obj){
  return(obj$fitted.values)
  
}

# coef method
coef.linreg <- function(obj){
  #names(obj$Coefficients)
  obj$Coefficients
}

# summary method
summary.linreg <- function(obj)
{
  se <- sqrt(diag(obj$vcov))
  tval <- obj$t_value
  estt <- obj$Coefficients
  TAB <- cbind(Estimate = estt,
               StdErr = se,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=obj$degfreedom))
  colnames(TAB) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
  res <- TAB
  #class(res) <- "summary.linreg"
  print(res)
  cat("\nResidual standard error:" , obj$varr ,"on", obj$degfreedom ,"degrees of freedom")
}


# plot method
library(cowplot)
plot.linreg <- function(object) {
  library(ggplot2)
  bp <- ggplot(data=iris,aes(x = object$fitted.values, y = object$residuals,label=Species))+ 
    geom_point(shape=1,size=4) + geom_smooth(method="lm", colour="red",
                                             se = FALSE)+ 
    labs(title="Residuals vs Fitted",x="Fitted values
      linreg(Petal.Length~Species, data = iris)", y="Residuals")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "white", colour = "grey50"))+
    
    geom_text(aes(label=ifelse(object$residuals>1&object$residuals<(-1),
                               as.character(Species),'')),hjust=0,vjust=0)
  
  sp <-   ggplot(data=iris,aes(x = object$fitted.values, y = object$standardized))+ 
    geom_point(shape=1,size=4) + geom_smooth(method="lm", colour="red",
                                             se = FALSE)+ 
    labs(title="Scale−Location",x="Fitted values
  linreg(Petal.Length~Species)", y=expression(sqrt(abs("Standard resduals"))))+
    theme(plot.title = element_text(hjust = 0.5),panel.background = 
            element_rect(fill = "white", colour = "grey50"))
  
  plot_grid(sp, bp, ncol = 1, nrow = 2)
  
}













