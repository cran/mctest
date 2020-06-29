\name{NEWS}
\title{News for Package 'mctest'}
\section{change in mctest version 1.3.1 (2020-06-29)}{
\itemize{
\item {Input argument} \code{x} in all functions of \code{mctest} package is changed to \code{model.matrix()} instead of \code{as.matrix()} for the acdeptance of factor variables as regressors. The numeric condition of \code{X} matrix is also omitted.
}
}


\section{Changes in mctest version 1.3 (2020-05-20)}{
\itemize{
\item {Input argument \code{x} and \code{y} are deferred.}
\item {A model object (\code{mod}), not necessarily type \code{lm} is taken as input argument in all of the functions \code{mctest()}, \code{imcdiag()}, \code{omcdiag()}, \code{eigprop()}, and \code{mc.plot()}.}
\item {VIF and Eigen values reference line (agains threshold value) in \code{mc.plot()} is enhanced.}
}
}


\section{Changes in mctest version 1.25 (2019-11-09)}{
\itemize{
\item {A Correction in is made in Farrar and Glauber Chi-Square formula. }
\item {Collinearity diagnostics \code{IND1} and \code{IND2} are added.}

}
}

\section{Changes in mctest version 1.2 (2018-02-24)}{
\itemize{
\item {A function \code{eigprop()} is added that computes the eigenvalues, condition indices and variance decomposition proportion collinearity diagnostic procedures found in Belsley, Kuh, and Welsch (1980).}
\item {References are updated and doi are added.}
\item {The x-axis label is removed to avoid the confusion of relating regressors with eigenvalues dimensions.}
}
}

\section{Changes in mctest version 1.11 (2017-10-28)}{
  \itemize{
  \item {References are updated.}
  \item {In \code{omcdiag()} and \code{imcdiag()} output, "COLLINEARITY in not detected by the test" is replaced by "COLLINEARITY is not detected by the test"}
  \item {In \code{omcdiag} output, "Condition Indices" is updated to "Condition Indices"}
}
}

\section{Changes in mctest version 1.1 (2016-12-05)}{
  \subsection{Addition}{
    \itemize{
      \item An option (\code{all}) for \code{imcdiag()} and \code{mctest()} is included that will show the output for all individual multicollinearity diagnostics as 1 or 0 if the values are above or below the threshold.
      \item Examples related to the use of \code{all} argument are added.
      \item Help files are improved.
    }
  }
  \subsection{Bug Fix}{
  \itemize{
    \item Bug is fixed.
  }
}
}
