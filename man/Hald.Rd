\name{Hald}
\alias{Hald}
\docType{data}
\title{Portland Cement benchmark of Hald(1952)}
\description{Heat evolved during setting of 13 cement mixtures of four basic ingredients. Each ingredient percentage appears to be rounded down to a full integer.  The sum of the four mixture percentages varies from a maximum of 99\% to a minimum of 95\%.  If all four regressor X-variables always summed to 100\%, the centered X-matrix would then be of rank only 3. Thus, the regression of heat on four X-percentages is ill-conditioned, with an approximate rank deficiency of MCAL = 1. The first column is the response and the remaining four columns are the predictors.

The Hald data as used by Hoerl, Kennard and Baldwin (1975). These data are also in package \code{wle}.
}
\usage{data(Hald)}
\format{
	A data frame with 13 observations on the following 5 variables.
	  \describe{
    \item{\code{Y}}{Heat (cals/gm) evolved in setting, recorded to nearest tenth.}
    \item{\code{X1}}{Integer percentage of 3CaO.Al2O3 in the mixture.}
    \item{\code{X2}}{Integer percentage of 3CaO.SiO2 in the mixture.}
    \item{\code{X3}}{Integer percentage of 4CaO.Al2O3.Fe2O3 in the mixture.}
    \item{\code{X4}}{Integer percentage of 2CaO.SiO2 in the mixture.}
  }

}

\source{
Woods H, Steinour HH, Starke HR. "Effect of composition of Portland cement on heat
 evolved during hardening. Industrial Engineering and Chemistry 1932; 24: 1207-1214.
}

\references{
Ridge Regression: some simulations, Hoerl, A. E. \emph{et al}, 1975, Comm Stat Theor Method 4:105
}
\examples{
data(Hald)
y <- Hald[ , 1]
x <- Hald[ , -1]
}
\keyword{datasets}
