\name{omcdiag}
\alias{omcdiag}
%\alias{print.omc}
%\alias{print.imc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Overall multicollinearity diagnostics}
\description{
Computes different overall measures of multicollinearit diagnostics for X matrix of regressors. Overall measures of collinearity detection includes Determinent of the correlation matrix (det), Farrar test of chi-square for presence of multicollinearity (Fchi), Red Indicator (Red), Sum of lambda inverse values (slambda), Theil's indicator and condition number with or without intercept term.
}
\usage{omcdiag(x, y, na.rm = TRUE, Inter = TRUE, detr=0.01, red=0.5,
                     conf=0.95, theil=0.5,cn=30,...) }

\arguments{
  \item{x}{a numeric design matrix and should contain more than one regressor}
  \item{y}{a numeric vector of response variable}
  \item{na.rm}{whether to remove missing observations}
  \item{Inter}{whether to include or exclude Intercept term}
  \item{detr}{determinant default threshold}
  \item{red}{red indicator default threshold}
  \item{conf}{confidence level of Farrar Chi-Square test}
  \item{theil}{Theil indicator default threshold}
  \item{cn}{condition number default threshold}
  \item{\dots}{extra argument(s) if used will be ignored}
}
\details{
This function detectes the existance of multicollinearity, by using different available diagnostic measure availble in literature such as Determinent of correlation matrix, Farrar test of chi-square, Red Indicator, Sum of lambda inverse values, Theil's Indicator and Condition Number.

Function also displays diagnostic measures value with the decision of either multicollinearity is detected by the diagnostics or not. The Value of 1 indicate that multicollineairty is detected and 0 indicate measure could not detect by the certain diagnostic measure.
}
\value{
  \item{idiags}{listing of overall diagnostic measures.  If \code{method} is not used all individual diagnostics will be dispalyed.}
  \item{Inter}{logical, if \code{TRUE} (the default value) a eigenvalues are returned with intercept term included.}
  \item{x}{matrix of regressors.}
}
\references{
Belsely, D.A., \emph{A Guide to Using the Collinearity Diagnostics}. Computer science in Economics and Management, 1991. \bold{4(1)}: 33--50.

Catterjee, S. and A. Hadi, \emph{Regression Analysis by Example}. 4 ed. 2006, Hoboken, New York; John Willey and Sons Inc.

Belsley, David. A., Edwin. Kuh, and Roy. E. Welsch. 1980. \emph{Regression Diagnostics: Identifying Influential Data and Sources of Collinearity}. New York: John Wiley and Sons.

Greene, William H. 2000. \emph{Econometric Analysis}. 4th Ed. Upper Saddle River, NJ: Prentice--Hall.

Imdadullah, M.  \emph{Addressing Linear Regression Models with Correlated Regressors: Some Package Development in R}, Unpublished Ph.D thesis, Department of Statistics, Bahauddin Zakriya University, Multan, Pakistan
}
\author{Muhammad Imdadullah, Muhammad Aslam}
\note{
Missing values in data will be removed by default. There is no method for the detection of multicollinearity, if missing values exists in the data set}

\seealso{
  \code{\link{imcdiag}}
}

\examples{
## Hald Cement data
data(Hald)
x<-Hald[,-1]
y<-Hald[,1]

## all oveall diagnostic measures and eigenvalues with intercept
od<-omcdiag(x,y)

## all oveall diagnostic measures and eigenvalues without intercept
omcdiag(x,y, Inter=FALSE)

## all oveall diagnostic measures and eigenvalues with intercept
## with different determinant and confidence level threshold

omcdiag(x,y, detr=0.001, conf=0.99)

## returns the determinant of correlation matrix |X'X|
omcdiag(x,y)[1]
}

\keyword{Farrar Chi-Square }
\keyword{Red Indicator}
\keyword{Condition Number}
\keyword{Theil's Indicator}
\keyword{Determinent of correlation matrix}
\keyword{sum of lambda inverse values}
