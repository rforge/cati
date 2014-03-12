\name{plot.ses.Tstats}
\alias{plot.ses.Tstats}

\title{
	Plot standardized effect size of T-statistics
}

\description{
	Plot standardized effect size of T-statistics
}

\usage{
	plot.ses.Tstats(tstats = NULL, val.quant = c(0.025, 0.975), col.Tstats = c("red", "purple", "green"), type = "normal", add.conf = TRUE)
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
	
	\item{type}{
		Type of plot. Possible type = "color_cond", "simple", "simple_sd", "normal" and "barplot".
	}
	
	\item{add.conf}{
		Logical value; Add confidence intervals or not. 
	}
}

\details{
	Normal type plot means, standard deviations, ranges and confidence intervals of T-statistics.
	Color_cond type plot T-statistics for each site and traits with color for significative values and grey for non signifivative ones. 
	Simple_sd type plot means, standard deviations and confidence intervals of T-statistics
	Simple type plot T-statistics for each site and traits and the mean confidence intervals by traits
	Barplot type plot means, standard deviations and confidence intervals of T-statistics in a barplot fashion
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
	\code{\link{plot.bar.Tstats}}
}

\examples{
	data(finch.ind)
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot_finch, sp=sp.finch, nperm=99)

	plot.ses.Tstats(res.finch)
	plot.ses.Tstats(res.finch, type="color_cond")
	plot.ses.Tstats(res.finch, type="simple")
	plot.ses.Tstats(res.finch, type="simple_sd")
	plot.ses.Tstats(res.finch, type="barplot")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
