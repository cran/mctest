\name{eigprop}
\alias{eigprop}
\alias{print.eigp}
\title{Eigenvalues and Variance Decomposition Proportion}

\description{
Computes eigenvalues, condition indices and variance decomposition proportions of \eqn{X'X}{X'X} or its related correlation matrix \eqn{R}{R} (see Belsley et al. (1980) <doi: 10.1007/BF00426854> ; Belsley, 1991; Kendall, 1957 and Silvey , 1969).
}

\usage{eigprop(x, na.rm = TRUE, Inter = TRUE, prop = 0.5, ...)}

\arguments{
  \item{x }{A numeric design matrix and should contain more than one regressor.}
  \item{na.rm }{Whether to remove missing observations.}
  \item{Inter}{Whether to include or exclude Intercept term, by default \code{Inter  =FALSE}.}
  \item{prop}{variance proportion default threshold, \code{prop=0.5}.}
  \item{\dots}{Extra argument(s) if used will be ignored.}
}
\details{
The \code{eigprop} function can be used to detect the existence of multicollinearity among regressors. The function computes eigenvalues, condition indices and variance decomposition proportions of regression coefficients. To check the linear dependencies associated with the corresponding eigenvalue, the \code{eigprop} compares variance proportion with threshold value (default is 0.5) and displays the proportions greater than given threshold from each row and column, if any. If \code{Inter = TRUE}, eigenvalues, condition indices and variance proportions are computed without intercept term. A list object of class "eigp" is returned:
}
\value{The \code{eigprop} objects are:
\item{ev}{A vector of eigenvalues. By default \code{Inter = TRUE}, eigenvalues are returned with intercept term included in the \code{X} matrix.}
\item{ci}{A vector of condition indices. By default \code{Inter = TRUE}, condition indices are returned with intercept term included in the \code{X} matrix.}
\item{call}{The matched call.}
\item{Inter}{logical, if \code{TRUE} (the default value) eigenvalues, condition indices and variance proportions are returned with intercept term included.}
\item{pi}{A matrix of variance decomposition proportions. By default \code{Inter = TRUE}, variance decomposition proportions are returned with intercept term included in the \code{X} matrix.}
\item{prop}{Default threshold proportion for comparison purpose.}
}

\references{
Belsely, D. A. A Guide to Using the Collinearity Diagnostics. \emph{Computer Science in Economics and Management}, \bold{4(1)}: 33--50, 1991.

Belsley, D. A., Kuh, E., and Welsch, R. E. \emph{Regression Diagnostics: Identifying Influential Data and Sources of Collinearity}. John Wiley \& Sons, New York, 1980.

Imdad, M. U. \emph{Addressing Linear Regression Models with Correlated Regressors: Some Package Development in R} (Doctoral Thesis, Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan), 2017.

Imdadullah, M., Aslam, M., and Altaf, S. mctest: An R Package for Detection of Collinearity Among Regressors. \emph{The R Journal}, \bold{8(2)}:499--509, 2016.

Silvey, S. D. Multicollinearity and imprecise estimation. \emph{Journal of the Royal Statistical Society, Series B (Methodological)}, \bold{31(3)}:539--552, 1969.
}
\author{Muhammad Imdad Ullah, Muhammad Aslam}

\note{Missing values in data will be removed by default. There is no method for the detection of multicollinearity, if missing values exists in the data set.}

\seealso{Overall collinearity diagnostics \code{\link{omcdiag}}, Individual collinearity diagnostics \code{\link{imcdiag}}}

\examples{
## Hald Cement data
data(Hald)
x<-Hald[,-1]

# with Intercept term
eigprop(x)

# without Intercept term
eigprop(x, Inter = FALSE)

# different proportion threshold
eigprop(x, prop = 0.45)

# only variance proportions
eigprop(x)$pi

# only condition indices
eigprop(x)$ci

# only eigenvalues
eigprop(x)$ev
}
\keyword{ eigenvalues }
\keyword{ condition number }
\keyword{ condition indices}
\keyword{ variance decomposition proportion}
\keyword{ multicollinearity }
