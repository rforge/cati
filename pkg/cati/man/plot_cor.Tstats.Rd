\name{plot_cor.Tstats}
\alias{plot_cor.Tstats}

\title{
	Plot relation between T-statistics
}

\description{
	Plot relation between the three T-statistics namely T_IP/IC, T_IC.IR and T_PC.PR.
}

\usage{
	plot_cor.Tstats(tstats = NULL, val.quant = c(0.025, 0.975), add.text= FALSE, bysites= FALSE, col.obj = NULL, plot.ask = TRUE, multipanel = TRUE, ...)
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
	
	\item{multipanel}{
		Logical value. If TRUE divides the device to shown several traits graphics in the same device. 
	}
	
	\item{\dots}{
		Any additional arguments are passed to the plot function creating the core of the plot and can be used to adjust the look of resulting graph.
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



\seealso{
	\code{\link{Tstats}};
	\code{\link{plot.Tstats}};
	\code{\link{plot_ses.var}}
}

\examples{
	data(finch.ind)
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot.finch, sp=sp.finch, nperm=9)
	
	plot_cor.Tstats(res.finch, bysites = FALSE)
	plot_cor.Tstats(res.finch, bysites = TRUE)
}




\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
