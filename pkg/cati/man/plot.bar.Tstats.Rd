\name{plot.bar.Tstats}
\alias{plot.bar.Tstats}

\title{
	Plot bar of T-statistics values
}

\description{
	Plot bar of T-statistics values.
}

\usage{
	plot.bar.Tstats(tstats = NULL, val.quant = c(0.025, 0.975), col.Tstats = c("red", "purple", "green", "white"))
}

\arguments{
	\item{tstats}{
		The list resulting from the function Tstats.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calcul confidence interval. By default val.quant=c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{col.Tstats}{
		A vector of three color correspond to the three T-statistics.
	}
}

\details{
}

\value{
None; used for the side-effect of producing a plot.
}

\references{
}

\author{
	Adrien Taudiere
}

\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{summary.Tstats}};
	\code{\link{Tstats}};
	\code{\link{plot.filter}};
	\code{\link{plot.ses}};
	\code{\link{plot.ses.Tstats}}
}

\examples{
	data(finch.ind)
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot_finch, sp=sp.finch, nperm=99)
	plot.bar.Tstats(res.finch)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
