\name{imcdiag}
\alias{imcdiag}
%\alias{print.omc}
%\alias{print.imc}
\title{Individual multicollinearity diagnostics}
\description{
Computes different measures of multicollinearity diagnostics for each regression in the design matrix X. Individual measures includes variance Inflation factor (VIF), Farrar F-test for determination of multicollinearity, Auxiliary F-test for relationship of F and R-square, Leamer's method, Corrected VIF (CVIF), and Klein's rule.

}
\usage{
imcdiag(x, y, method = NULL, na.rm = TRUE, corr = FALSE, vif=10, tol=0.1,
              conf=0.95, cvif=10, leamer=0.1,...)
}
\arguments{
  \item{x }{a numeric design matrix and should contain more than one regressor}
  \item{y }{a numeric vector of response variable}
  \item{na.rm }{whether to remove missing observations}
  \item{method }{specific individual measure of collinearity such as VIF, CVIF, and Leamer, etc. For example, \code{method="VIF"}}
  \item{corr }{whether to display correlation matrix or not \code{Inter=TRUE}}
  \item{vif }{defualt threshold for VIF measure, \code{vif=10}}
  \item{tol}{default thereshold for TOL measure, \code{tol=0.10}}
  \item{conf }{default confidence level for Farrar's Wi test, \code{conf=0.99}}
  \item{cvif }{default threshold for CVIF measure, \code{CVIF=10}}
  \item{leamer }{default threshold for Leamer's method, \code{leamer=0.1}}
  \item{...}{extra argument(s) if used will be ignored}
}
\details{
The \code{imcdiag} function detectes the existance of multicollinearity due to x-variable. That's why named as individual measures of diagnostics. This includes VIF, TOL, Condition Indices, Klien's rule, Theil's measure, Red Indicator, Farrar and Glauber test of Chi-square and F-test, Leamer's method etc. If \code{method} argument is used (\code{method="VIF"}), the result for VIF indicator will be displayed with decision of either collinearity exists or not due to certain regressor which is indicated by 0 (collinearity is not detected by method for regressor given in row) and 1 (collinearity is detected by the method for regressor given in row).

}

\value{
Returns the diagnostic measures value with decision of either multicollinearity is detected by the diagnostic measure or not. Value of 1 indicate multicollineairty is detected and 0 indicate measure could not detect the existance of multicollinearity. An output list object of class imc:
  \item{idiags}{listing of specific individual measure such as \code{method="CVIF"} provided.  If \code{method} is not used all individual diagnostics will be dispalyed.}
  \item{x}{matrix of regressors.}
  \item{y}{a vector of response variable.}
  \item{method}{individual collinearity measure.}
  \item{corr}{logical, if \code{FALSE} (the default value) a correlation matrix is returned after all or required individual collinearity diagnostic.}
}
\references{
Belsely, D.A., \emph{A Guide to Using the Collinearity Diagnostics}. Computer science in Economics and Management, 1991. \bold{4(1)}: 33--50.

Catterjee, S. and A. Hadi, \emph{Regression Analysis by Example}. 4 ed. 2006, Hoboken, New York; John Willey and Sons Inc.

Belsley, David. A., Edwin. Kuh, and Roy. E. Welsch. 1980. \emph{Regression Diagnostics: Identifying Influential Data and Sources of Collinearity}. New York: John Wiley and Sons.

Greene, William H. 2000. \emph{Econometric Analysis}. 4th Ed. Upper Saddle River, NJ: Prentice--Hall.
}
\author{Muhammad Imdadullah, Dr. Muhammad Aslam}
\note{
Missing values in data will be removed by default. There is no method for the detection of multicollinearity, if missing values exists in the data set}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{omcdiag}} }
\examples{
## Hald Cement data
data(Hald)
x<-Hald[,-1]
y<-Hald[,1]

## all Individual measures
id<-imcdiag(x,y); id$idiags[,1]

# VIF measure with custom VIF threshold
imcdiag(x,y, method="VIF", vif=5)

# CVIF measure with custom CVIF threshold and correlation matrix
imcdiag(x,y, method="CVIF", cvif=5, corr=TRUE)

## only VIF values without collinearity detection indication
imcdiag(x,y, method="VIF")[[1]][,1]
plot(imcdiag(x,y, method="VIF")[[1]][,1]) # vif plot
}
\keyword{VIF}
\keyword{CVIF}
\keyword{TOL}
\keyword{Farrar F-test}
\keyword{F-test and R-Square relation}
\keyword{klein's Rule}
\keyword{Leamer's Method}
