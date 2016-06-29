\name{mctest}
\alias{mctest}

\title{ Multicollinearity diagnostic measures }
\description{
\code{mctest} display overall, individual or both types of multicollinearity diagnostic measures from      \code{omcdiag} and \code{imcdiag} functions respectively.
}
\usage{
mctest(x, y, type=c("o","i","b"), na.rm = TRUE, Inter=TRUE, method=NULL,
        corr=FALSE, detr=0.01, red=0.5, theil=0.5, cn=30, vif=10, tol=0.1,
        conf=0.95, cvif=10, leamer=0.1, ...)
}
\arguments{
  \item{x}{a numeric design matrix and should contain more than one regressor}
  \item{y}{a numeric vector of response variable}
  \item{na.rm}{whether to remove missing observations}
  \item{Inter}{whether to include or exclude Intercept term}
  \item{type}{displays overall, individual or both type of collinearity dignostics. Individual collinearity diagnostics are displayed by default.}
  \item{method }{specific individual measure of collinearity such as VIF, TOL, CVIF, adn Leamer, etc. For example, \code{method="VIF"}}
  \item{corr }{whether to display correlation matrix or not \code{Inter=TRUE}}
\item{detr}{determinant default threshold}
  \item{red}{red indicator default threshold}
  \item{theil}{Theil indicator default threshold}
  \item{cn}{condition number default threshold}
  \item{vif }{defualt threshold for VIF measure, \code{vif=10}}
  \item{conf }{default confidence level for Farrar's test, \code{conf=0.99}}
  \item{cvif }{default threshold for CVIF measure, \code{CVIF=10}}
  \item{tol }{default threshold forTOL measure, \code{TOL=0.10}}
  \item{leamer }{default threshold for Leamer's method, \code{leamer=0.1}}
  \item{\dots}{extra argument(s) if used will be ignored}
}

\references{
Imdadullah, M.  \emph{Addressing Linear Regression Models with Correlated Regressors: Some Package Development in R}, Unpublished Ph.D thesis, Department of Statistics, Bahauddin Zakriya University, Multan, Pakistan.

Belsely, D.A., \emph{A Guide to Using the Collinearity Diagnostics}. Computer science in Economics and Management, 1991. \bold{4(1)}: 33--50.

Catterjee, S. and A. Hadi, \emph{Regression Analysis by Example}. 4 ed. 2006, Hoboken, New York; John Willey and Sons Inc.

Belsley, David. A., Edwin. Kuh, and Roy. E. Welsch. 1980. \emph{Regression Diagnostics: Identifying Influential Data and Sources of Collinearity}. New York: John Wiley and Sons.

Greene, William H. 2000. \emph{Econometric Analysis}. 4th Ed. Upper Saddle River, NJ: Prentice--Hall.

}
\author{Muhammad Imdadullah, Dr. Muhammad Aslam}

\note{ Missing values in data will be removed by default. There is no method for the detection of multicollinearity, if missing values exists in the data set}

\seealso{\code{\link{omcdiag}, \link{imcdiag}}}
\examples{
## Hald Cement data
data(Hald)
x<-Hald[,-1]
y<-Hald[,1]

## Overall diagnostic measures and eigenvalues with intercept term
mctest(x, y)

## Overall diagnostic measures and eigenvalues without intercept term
mctest(x, y, Inter=FALSE)

## all individual diagnostic measures
mctest(x, y, type="i")

## certain individual diagnostic measures with collinearity detection indication
VIF<-mctest(x, y, type="i", method="VIF")
VIF[[1]][,1] # named VIF values only

## all individual diagnostic measures with correlation matrix
mctest(x, y, type="i", corr=TRUE)

## VIF and correlation matrix with collinearity detection indication
mctest(x, y, type="i", method="VIF", corr=TRUE)

## both overall and individual collinearity diagnostics
mctest(x, y, type="b")
mctest(x, y, type="b", method="VIF", cor=TRUE)

## all overall and vif with correlation matrix
## VIF and CN desired threshold
## eigenvalues without intercept term
mctest(x, y, type="b", method="VIF", Inter=FALSE, vif=15, cn=35)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
