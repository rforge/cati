\name{plot.cor.Tstats}
\alias{plot.cor.Tstats}

\title{
	Plot relation between T-statistics
}

\description{
	Plot relation between the three T-statistics namely T_IP/IC, T_IC.IR and T_PC.PR.
}

\usage{
plot.cor.Tstats(tstats = NULL, val.quant = c(0.025, 0.975), add.text = FALSE, bysites = FALSE, col.obj = NULL, plot.ask = TRUE)
}

\arguments{
	\item{tstats}{
		The list resulting from the function Tstats.
	}
	
	\item{val.quant}{		
		Numeric vectors of length 2, giving the quantile to calcul confidence interval. By default val.quant=c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{add.text}{
		Logical value; Add text or not.
	}
	
	\item{bysites}{
		Logical value; plot by sites or by traits. 
	}
	
	\item{col.obj}{
		Vector of colors for object (either traits or sites).
	}
	
	\item{plot.ask}{
		Logical value; Ask for new plot or not.
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
	\code{\link{Tstats}};
	\code{\link{plot.ses.Tstats}};
	\code{\link{plot.ses.var}}
}

\examples{
	data(finch.ind)
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot_finch, sp=sp.finch, nperm=99)
	
	plot.cor.Tstats(res.finch, bysites = FALSE)
	plot.cor.Tstats(res.finch, bysites = TRUE)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
