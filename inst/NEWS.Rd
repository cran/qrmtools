% Check from R:
% news(db = tools:::.build_news_db_from_package_NEWS_Rd("NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{qrmtools}}
\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in version 0.0-10}{
  \subsection{New Features}{
    \itemize{
      \item \code{dGPDtail()}, \code{pGPDtail()}, \code{qGPDtail()},
      \code{rGPDtail()} (and renamed \code{VaR_POT()} to
      \code{VaR_GPDtail()} and \code{ES_POT()} to \code{ES_GPDtail()}
      correspondingly)
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{plot_matrix(): used to plot transposed matrix.}
    }
  }
}

\section{Changes in version 0.0-9}{
  \subsection{New Features}{
    \itemize{
      \item \code{VaR_POT()}, \code{ES_POT()}
      \item \code{GPD_shape_plot()}, \code{tail_plot()}
      \item \code{mean_excess_np()}, \code{mean_excess_GPD()} and
            \code{mean_excess_plot()}
      \item \code{edf_plot()}
      \item Fitting of GEV distribution and GPD.
      \item More intuitive (re)naming of all arguments of EVT related functions.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item (none)
    }
  }
}
