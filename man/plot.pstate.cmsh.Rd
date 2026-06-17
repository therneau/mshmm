\name{plot.pstate.cmsh}
\alias{plot.pstate.cmsh}
\title{
Plot predicted probability in state curves
}
\description{
  This is the plot method for pstate.cmsh objects produced by the
  \code{predict.cmsh} method.
}
\usage{
\method{plot}{pstate.cmsh}(x, cumhaz = FALSE, cumprob = FALSE, 
    xlab, ylab, xscale = 1, yscale = 1, xmax, xlim, ylim, ...)
}
\arguments{
  \item{x}{an object of class \code{pstate.cmsh}, usually created by
    the \code{predict.cmsh} function
}
  \item{cumhaz}{plot the cumulative hazard rather than the probability
    in state or survival.  Optionally, this can be a numeric vector
    specifying which columns of the \code{cumhaz} component to plot.
}
  \item{cumprob}{for multi-state curves, plot the probabilities in
    state 1, (state1 + state2), (state1 + state2 + state3), \ldots.
    If \code{cumprob} is an integer vector the totals will be in the
    order indicated, a logical vector sums the selected states
}
  \item{xlab, ylab}{x and y axis labels}
  \item{xscale}{optional scaling for the horizontal axis. A common
    choice when the curve's timscale is days is \code{xcale=365.25}.
    This only changes the labels, not the underlying scale of the data;
    subsequent actions like \code{locator()} are on the original scale.
}
  \item{yscale}{
    a numeric value used to multiply the labels on the y axis. 
    A value of 100, for instance, would be used to give a percent scale. 
    Only the labels are 
    changed, not the actual plot coordinates, so that adding a curve
    from another fit 
    will perform as it did without the \code{yscale} argument. 
  }
  \item{xmax}{
    the maximum horizontal plot coordinate.  This can be used to shrink
    the range of a plot.  It shortens the curve before plotting it, so
    unlike using the \code{xlim} graphical parameter, warning
    messages about out of bounds points are not generated.
    It also \emph{looks} nicer when lines don't touch the right margin.
}
\item{xlim}{optional limit for the x axis}
\item{ylim}{optional limit for the  axis}
  \item{\dots}{other graphics arguments such as lty, col, lwd, etc.
}
}
\keyword{ graphics }
