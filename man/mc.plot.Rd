\name{mc.plot}
\alias{mc.plot}
\title{Plot of VIF and Eigen values}
\description{
Plot of VIF and Eigen values for detection of multicollinearity among regressors. The VIF and Eigen values are also displayed on graph. Eigen values plot can be displayed with or without inclusion of intercept term.
}
\usage{mc.plot(x,y, Inter=FALSE, vif=10, ev=0.01, ...)}
\arguments{
  \item{x}{a numeric design matrix and should contain more than one regressor}
  \item{y}{a numeric vector of response variable}
  \item{Inter}{whether to include or exclude Intercept term}
  \item{vif}{threshold of VIF and will appear as horizontal line on VIF plot. The default value is \code{vif=10}.}
  \item{ev}{threshold of Eigenvalues and will appear as horizontal line on Eigvenvalues plot. The default value is \code{ev=0.01}.}
  \item{\dots}{extra argument(s) if used will be ignored}
}
\details{
\code{mc.plot} function draw graphs of VIF and Eigenvalues for graphical detection of collinearity among regression. Horizontal line for VIF and Eigenvalues is drawn as indication of multicollinearity existence.
}
\value{Don't return any thing, it displays plot}

\references{Imdadullah, M.  \emph{Addressing Linear Regression Models with Correlated Regressors: Some Package Development in R}, Unpublished Ph.D thesis, Department of Statistics, Bahauddin Zakriya University, Multan, Pakistan.

Belsely, D.A., \emph{A Guide to Using the Collinearity Diagnostics}. Computer science in Economics and Management, 1991. \bold{4(1)}: 33--50.

Catterjee, S. and A. Hadi, \emph{Regression Analysis by Example}. 4 ed. 2006, Hoboken, New York; John Willey and Sons Inc.

Belsley, David. A., Edwin. Kuh, and Roy. E. Welsch. 1980. \emph{Regression Diagnostics: Identifying Influential Data and Sources of Collinearity}. New York: John Wiley and Sons.

Greene, William H. 2000. \emph{Econometric Analysis}. 4th Ed. Upper Saddle River, NJ: Prentice--Hall.

}

\author{
Muhammad Imdadullah, Muhammad Aslam
}

\seealso{\code{\link{omcdiag}, \link{imcdiag}}}

\examples{
## Hald Cement data
data(Hald)
x<-Hald[,-1]
y<-Hald[,1]

## plot with default threshold of VIF and Eigenvalues with no intercept
mc.plot(x, y)

## plot with default threshold of VIF and Eigenvalues with intercept
mc.plot(x, y, Inter=TRUE)

## plot with specific threshold of VIF and Eigenvalues with no intercept
mc.plot(x, y, vif=5, ev=20)

## plot with specific threshold of VIF and Eigenvalues with intercept
mc.plot(x, y, vif=5, ev=20, Inter=TRUE)
}
\keyword{VIF }
\keyword{Eigenvalues}
