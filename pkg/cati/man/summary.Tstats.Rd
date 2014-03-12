\name{summary.Tstats}
\alias{summary.Tstats}

\title{
	Summarize T-statistics
}

\description{
	Function to summarize t statistics
}

\usage{
	summary.Tstats(tstats = NULL, val.quant = c(0.025, 0.975), type = "all")
}

\arguments{
	\item{tstats}{
		The list resulting from the function Tstats.
	}

	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calcul confidence interval. By default val.quant=c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{type}{
		Type of summary statistics. Either 'binary', 'percent', 'p.value', 'site' or 'all'. 
	}
	
}

\details{
	Binary type only test if a T-statistics is significatively different from the null expectation for each traits.
	Percent type determine the percentage of sites were the T-statistics is significatively different from the null expectation for each traits. Asterix shows global significance of the test.
	P-value type determine the p-value (two unilateral tests) of the T-statistics for each traits and sites.
	Site type allows to know in which sites T-statistics deviate from the null expectation.
	All type do all the precedent type of summary.
}

\value{
	A matrix or a list of matrix if type=all.
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
	\code{\link{Tstats}};
	\code{\link{plot.ses.Tstats}};
	\code{\link{plot.filter}};
	\code{\link{mat.comm}}
}

\examples{
	data(finch.ind)
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot_finch, sp=sp.finch, nperm=99)
	
	summary.Tstats(res.finch, type="binary")
	summary.Tstats(res.finch, type="percent")
	summary.Tstats(res.finch, type="site")
	summary.Tstats(res.finch, type="p.value")
	summary.Tstats(res.finch, type="all")
	
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
