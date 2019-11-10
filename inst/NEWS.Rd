\name{NEWS}
\title{News for Package 'mctest'}
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
