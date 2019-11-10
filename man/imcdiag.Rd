\name{imcdiag}
\alias{imcdiag}
%\alias{print.omc}
\alias{print.imc}
\title{Individual Multicollinearity Diagnostic Measures}
\description{
Computes different measures of multicollinearity diagnostics for each regressor in the design matrix \eqn{X}{X}. Individual measures includes variance Inflation factor (VIF) (Marquardt, 1970), Farrar F-test for determination of multicollinearity (Farrar and Glauber, 1967), Auxiliary F-test for relationship between F and R-square(Gujarati and Porter, 2008), Leamer's method (Greene, 2002), Corrected VIF (CVIF) Curto and Pinto (2011) <doi: 10.1080/02664763.2010.505956>, Klein's rule Klein (1962), and IND1 & IND2 (Imdad, et. al., 2019) <\url{https://doi.org/10.17576/jsm-2019-4809-26}> proposed by the researchers.

}
\usage{
imcdiag(x, y, method = NULL, na.rm = TRUE, corr = FALSE, vif = 10, tol = 0.1,
              conf = 0.95, cvif = 10, ind1 = 0.02, ind2 = 0.7, leamer = 0.1,
              all = FALSE, ...)
}
\arguments{
  \item{x }{A numeric design matrix and should contain more than one regressor.}
  \item{y }{A numeric vector of response variable.}
  \item{na.rm }{Whether to remove missing observations.}
  \item{method }{Specific individual measure of collinearity such as VIF, CVIF, and Leamer, etc. For example, \code{method="VIF"}.}
  \item{corr }{Whether to display correlation matrix or not, by default \code{corr=FALSE}.}
  \item{vif }{Default threshold for VIF measure, \code{vif=10}.}
  \item{tol}{Default threshold for TOL measure, \code{tol=0.10}.}
  \item{conf }{Default confidence level for Farrar's Wi test, \code{conf=0.99}.}
  \item{cvif }{Default threshold for CVIF measure, \code{CVIF=10}.}
  \item{ind1 }{Default threshold for IND1 indicator, \code{ind1=0.02}}
  \item{ind2 }{Default threshold for IND2 indicator, \code{ind2=0.7}}
  \item{leamer }{Default threshold for Leamer's method, \code{leamer=0.1}.}
  \item{all }{Returns all individual measure of collinearity in a matrix of 0 (not detected) or 1 (detected).}
  \item{...}{Extra argument(s) if used will be ignored.}
}
\details{
The \code{imcdiag} function detects the existence of multicollinearity due to \eqn{x}{x}-variable. That's why named as individual measures of diagnostics. This includes VIF, TOL, Klein's rule, Farrar and Glauber F-test, F and \eqn{R^2}{R^2} relation, Leamer's method, CVIF, IND1, and IND2 diagnostic measures of multicollinearity. If \code{method} argument is used (\code{method="VIF"}), the VIF values for each regressor will be displayed with decision of either collinearity exists or not which is indicated by 0 (collinearity is not detected by method for regressor) and 1 (collinearity is detected by the method for regressor). If argument \code{all=TRUE} all individual measures of collinearity will be displayed in a matrix of 0 (collinearity is not detected) or 1 (collinearity is detected).
}

\value{
This function detects the existence of multicollinearity by using different available diagnostic measures already available in literature. The function returns the value of diagnostic measures with decision of either collinearity is detected by the diagnostic measure or not. Value of 1 indicates that collinearity is detected and 0 indicates that measure could not detect the existence of collinearity. A list object of class "imc" is returned:
  \item{idiags}{Listing of specific individual measure such as \code{method="CVIF"} provided. If \code{method} is not used all individual diagnostics will be displayed.}
  \item{x}{A numeric matrix of regressors.}
  \item{y}{A vector of response variable.}
  \item{method}{Specific individual collinearity measure, such as VIF, TOL, CVIF, IND1, and IND2 etc.}
  \item{corr}{Logical, if \code{FALSE} (the default value) a correlation matrix will not be displayed.}
  \item{R2}{R-square from regression of all regressors \eqn{X}{X} on response variable \eqn{y}{y}.}
  \item{call}{The matched call.}
  \item{pval}{Returns significant regressor as number after comparing the \emph{p}-value of regressors from \code{summary.lm} function with \eqn{1-conf}{1-conf}. }
  \item{all}{If \code{TRUE} individual collinearity measures will be returned as a matrix of 0 or 1.}
  \item{alldiag}{Matrix of all individual collinearity measures indicated as either 0 (collinearity not detected) or 1 (collinearity detected) for each diagnostic measure and each regressor.}
}

\references{
Belsely, D. A. A Guide to Using the Collinearity Diagnostics. \emph{Computer Science in Economics and Management}, \bold{4(1)}: 33--50, 1991.

Belsley, D. A., Kuh, E., and Welsch, R. E. \emph{Regression Diagnostics: Identifying Influential Data and Sources of Collinearity}. John Wiley \& Sons, New York, 1980.

Chatterjee, S. and Hadi, A. S. \emph{Regression Analysis by Example}. John Wiley \& Sons, 4th edition, New York, 2006.

Curto, J. D. and Pinto, J. C. The Corrected VIF (CVIF). \emph{Journal of Applied Statistics}, \bold{38}(7), 1499--1507.

Greene, W. H. \emph{Econometric Analysis}. Prentice--Hall, Upper Saddle River, New Jersey, 4th edition, 2000.

Imdad, M. U. \emph{Addressing Linear Regression Models with Correlated Regressors: Some Package Development in R} (Doctoral Thesis, Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan), 2017.

Imdadullah, M., Aslam, M., and Altaf, S. mctest: An R Package for Detection of Collinearity Among Regressors. \emph{The R Journal}, \bold{8}(2):499--509, 2016.

Imdad, M. U., Aslam, M., Altaf, S., and Ahmed, M. Some New Diagnostics of Multicollinearity in Linear Regression Model. \emph{Sains Malaysiana}, \bold{48}(2): 2051--2060, 2019.
}

\author{Muhammad Imdad Ullah, Muhammad Aslam}

\note{Missing values in data will be removed by default. There is no method for the detection of multicollinearity, if missing values exists in the data set.}

\seealso{Overall collinearity diagnostic \code{\link{omcdiag}}, collinearity plot \code{\link{mc.plot}} }

\examples{
## Hald Cement data
data(Hald)
x<-Hald[,-1]
y<-Hald[,1]

## all Individual measures
id<-imcdiag(x, y); id$idiags[,1]

# VIF measure with custom VIF threshold
imcdiag(x, y, method = "VIF", vif = 5)

# IND1 measure with custom IND1 threshold and correlation matrix
imcdiag(x,y, method="IND1", ind1=0.01, corr=TRUE)

# CVIF measure with custom CVIF threshold and correlation matrix
imcdiag(x, y, method = "CVIF", cvif = 5, corr = TRUE)

# Collinearity Diagnostic measures in matrix of 0 or 1
imcdiag(x, y, all = TRUE)
imcdiag(x, y, method = "VIF", all = TRUE)

## only VIF values without collinearity detection indication
imcdiag(x, y, method = "VIF")[[1]][,1]
plot(imcdiag(x, y, method = "VIF")[[1]][,1]) # vif plot
}
\keyword{VIF}
\keyword{CVIF}
\keyword{TOL}
\keyword{Farrar F-test}
\keyword{F-test and R-Square relation}
\keyword{IND1}
\keyword{IND2}
\keyword{klein's Rule}
\keyword{Leamer's Method}
\keyword{Multicollinearity}
\keyword{Leamer's Method}
