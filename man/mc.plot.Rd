\name{mc.plot}
\alias{mc.plot}
\title{Plot of VIF and Eigenvalues}
\description{
Plot of VIF and Eigenvalues for detection of multicollinearity among regressors. The VIF and Eigenvalues are also displayed on graph. Eigenvalues plot can be displayed with or without inclusion of intercept term.
}
\usage{mc.plot(x, y, Inter = FALSE, vif = 10, ev = 0.01, ...)}
\arguments{
  \item{x}{A numeric design matrix and should contain more than one regressor.}
  \item{y}{A numeric vector of response variable.}
  \item{Inter}{Whether to include or exclude Intercept term, by default \code{Inter=FALSE}.}
  \item{vif}{Threshold of VIF and will appear as horizontal line on VIF plot. The default value is \code{vif=10}.}
  \item{ev}{Threshold of Eigenvalues and will appear as horizontal line on Eigenvalues plot. The default value is \code{ev=0.01}.}
  \item{\dots}{Extra argument(s) if used will be ignored.}
}
\details{
\code{mc.plot} function draw graphs of VIF and Eigenvalues for graphical detection of collinearity among regression. Horizontal line for VIF and Eigenvalues is drawn as indication of threshold values of both VIF and Eigenvalues for testing the multicollinearity.
}
\value{Don't return any thing, it displays plot.}

\references{
Belsely, D. A. A Guide to Using the Collinearity Diagnostics. \emph{Computer Science in Economics and Management}, \bold{4(1)}: 33--50, 1991.

Belsley, D. A., Kuh, E., and Welsch, R. E. \emph{Regression Diagnostics: Identifying Influential Data and Sources of Collinearity}. John Wiley \& Sons, New York, 1980.

Chatterjee, S. and Hadi, A. S. \emph{Regression Analysis by Example}. John Wiley \& Sons, 4th edition, New York, 2006.

Greene, W. H. \emph{Econometric Analysis}. Prentice--Hall, Upper Saddle River, New Jersey, 4th edition, 2000.

Imdad, M. U. \emph{Addressing Linear Regression Models with Correlated Regressors: Some Package Development in R} (Doctoral Thesis, Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan), 2017.

Imdadullah, M., Aslam, M., and Altaf, S. mctest: An R Package for Detection of Collinearity Among Regressors. \emph{The R Journal}, \bold{8(2)}:499--509, 2016.
}

\author{Dr. Muhammad Imdad Ullah, Dr. Muhammad Aslam}

\seealso{Overall collinearity diagnostic measures \code{\link{omcdiag}, Individual collinearity diagnostic measures \link{imcdiag}}}

\examples{
## Hald Cement data
data(Hald)
x<-Hald[,-1]
y<-Hald[,1]

## plot with default threshold of VIF and Eigenvalues with no intercept
mc.plot(x, y)

## plot with default threshold of VIF and Eigenvalues with intercept
mc.plot(x, y, Inter = TRUE)

## plot with specific threshold of VIF and Eigenvalues with no intercept
mc.plot(x, y, vif = 5, ev = 20)

## plot with specific threshold of VIF and Eigenvalues with intercept
mc.plot(x, y, vif = 5, ev = 20, Inter = TRUE)
}
\keyword{VIF }
\keyword{Eigenvalues}
