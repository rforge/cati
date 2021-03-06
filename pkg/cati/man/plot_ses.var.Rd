\name{plot_ses.var}
\alias{plot_ses.var}

\title{
	Plot SES values against a variable
}

\description{
	Plot standardized effect size values against a variable
}

\usage{
	plot_ses.var(index.list, variable = NULL, ylab = "variable", color.traits = NULL, val.quant = c(0.025, 0.975), resume = FALSE, multipanel = TRUE)
}

\arguments{
	\item{index.list}{
		A list of index and the associate null models in the forme: list( index_1=index_1_observed, index_1_nm=null.model.index_1 ,index_2=index_2_observed, index_2_nm=null.model.index_2, ...). 
	}
	
	\item{variable}{
		The variable to plot standardized effect size of index.
	}
	
	\item{ylab}{
		Label for the variable.
	}
	
	\item{color.traits}{
		A vector of colors corresponding to traits.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calcul confidence interval. By default val.quant=c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{resume}{
		Logical value; resume=FALSE by default; Simplify the plot by ploting the mean and standard error for index value of multiple traits
	}
	
	\item{multipanel}{
		Logical value. If TRUE divides the device to shown several traits graphics in the same device. 
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
	\code{\link{plot.listofindex}};
	\code{\link{ses}}
}

\examples{
	data(finch.ind)
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot.finch, sp=sp.finch, nperm=9)

	par(mfrow=c(2,2))
	species.richness<-table(ind.plot.finch)
	plot_ses.var(as.listofindex(list(res.finch)), species.richness, multipanel=FALSE)

	#Same plot with \code(resume=TRUE).
	
	par(mfrow=c(2,2))
	plot_ses.var(as.listofindex(list(res.finch)), species.richness, resume=TRUE, multipanel=FALSE)
	par(mfrow=c(1,1))
}


\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
