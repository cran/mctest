\name{mctest}
\alias{mctest}

\title{ Multicollinearity diagnostic measures }

\description{
The function \code{mctest} display overall, individual or both types of multicollinearity diagnostic measures from \code{omcdiag} and \code{imcdiag} functions, respectively.}

\usage{
mctest(x, y, type=c("o","i","b"), na.rm = TRUE, Inter=TRUE, method=NULL,
        corr=FALSE, detr=0.01, red=0.5, theil=0.5, cn=30, vif=10, tol=0.1,
        conf=0.95, cvif=10, ind1=0.02, ind2=0.7, leamer=0.1, all=FALSE, ...)
}
\arguments{
  \item{x}{A numeric design matrix and should contain more than one regressor.}
  \item{y}{A numeric vector of response variable.}
  \item{na.rm}{Whether to remove missing observations.}
  \item{Inter}{Whether to include or exclude Intercept term. By default \code{Inter=TRUE}.}
  \item{type}{Displays overall, individual or both type of collinearity diagnostics. Overall collinearity diagnostics are displayed by default with eigenvalues and condition indexes, when \code{method} and \code{type} argument are not used.}
  \item{method }{Specific individual measure of collinearity such as VIF, TOL, CVIF, Leamer, IND1, and IND2  etc, when \code{method} argument is used. For example, \code{method="VIF"}.}
  \item{corr }{Whether to display correlation matrix or not \code{Inter=TRUE}.}
\item{detr}{Determinant default threshold, \code{detr=0.01}.}
  \item{red}{Red indicator default threshold, \code{red=0.5}.}
  \item{theil}{Theil's indicator default threshold, \code{theil=0.5}.}
  \item{cn}{Condition number default threshold, \code{cn=30}.}
  \item{vif }{Default threshold for VIF measure, \code{vif=10}.}
  \item{conf }{Default confidence level for Farrar's test, \code{conf=0.99}.}
  \item{cvif }{Default threshold for CVIF measure, \code{CVIF=10}.}
  \item{tol }{Default threshold for TOL measure, \code{TOL=0.10}.}
  \item{ind1 }{Default threshold for IND1 indicator, \code{ind1=0.02}.}
  \item{ind2 }{Default threshold for IND2 indicator, \code{ind2=0.7}.}
  \item{leamer }{Default threshold for Leamer's method, \code{leamer=0.1}.}
  \item{all }{Returns all individual measure of collinearity in a matrix of 0 (not detected) or 1 (detected).}
  \item{\dots}{Extra argument(s) if used will be ignored.}
}

\references{
Belsely, D. A. A Guide to Using the Collinearity Diagnostics. \emph{Computer Science in Economics and Management}, \bold{4(1)}: 33--50, 1991.

Belsley, D. A., Kuh, E., and Welsch, R. E. \emph{Regression Diagnostics: Identifying Influential Data and Sources of Collinearity}. John Wiley \& Sons, New York, 1980.

Chatterjee, S. and Hadi, A. S. \emph{Regression Analysis by Example}. John Wiley \& Sons, 4th edition, New York, 2006.

Greene, W. H. \emph{Econometric Analysis}. Prentice--Hall, Upper Saddle River, New Jersey, 4th edition, 2000.

Imdad, M. U. \emph{Addressing Linear Regression Models with Correlated Regressors: Some Package Development in R} (Doctoral Thesis, Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan), 2017.

Imdadullah, M., Aslam, M., and Altaf, S. mctest: An R Package for Detection of Collinearity Among Regressors. \emph{The R Journal}, \bold{8(2)}:499--509, 2016.

Imdad, M. U., Aslam, M., Altaf, S., and Ahmed, M. Some New Diagnostics of Multicollinearity in Linear Regression Model. \emph{Sains Malaysiana}, \bold{48}(2): 2051--2060, 2019.
}

\author{Muhammad Imdad Ullah, Muhammad Aslam}

\note{ Missing values in data will be removed by default. There is no method for the detection of multicollinearity, if missing values exists in the data set}

\seealso{overall collinearity diagnostics \code{\link{omcdiag}, individual collinearity diagnostics \link{imcdiag}}, collinearity plots \code{\link{mc.plot}}}

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

IND1<-mctest(x, y, type="i", method="IND1")
IND1[[1]][,1] # named IND1 values only

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

## Individual collinearity diagnostic measures in matrix of 0 or 1
mctest(x, y, all = TRUE)
mctest(x, y, method = "VIF", all = TRUE)
mctest(x, y, type="b", all = TRUE)

}

\keyword{ Multicollinearity }
\keyword{ Overall collinearity diagnostic measures }
\keyword{ Individual collinearity diagnostic measures }
