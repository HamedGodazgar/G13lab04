\name{linreg}
\alias{linreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
linreg(formula, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, data) 
{
    X <- model.matrix(formula, data)
    y <- data[all.vars(formula)[1]]
    QR <- qr(X)
    r <- QR$rank
    piv <- QR$pivot[1:r]
    betahat <- solve.qr(QR, y)
    betahat
    yhat <- as.vector(c(X[, piv] \%*\% betahat))
    y <- as.matrix(data[all.vars(formula)[1]])
    resi <- as.vector(y - yhat)
    dof <- nrow(X) - r
    se2 <- c(crossprod(resi))/dof
    v <- chol2inv(QR$qr, r) * se2
    tval <- betahat/sqrt(diag(v))
    yy <- data[all.vars(formula)[1]]
    return_obj <- list(Coefficients = betahat, fitted.values = yhat, 
        residuals = resi, degfreedom = dof, vcov = v, varr = sqrt(se2), 
        t_value = tval, call = match.call(), y = yy)
    class(return_obj) <- "linreg"
    return(return_obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
