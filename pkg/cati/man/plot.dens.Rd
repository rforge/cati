\name{plot.dens}
\alias{plot.dens}

\title{
	Plot function to represent density of traits values
}

\description{
	Plot function to represent density of traits values
}

\usage{
	plot.dens(traits = NULL, var.1 = NULL, var.2 = NULL, col.dens = NULL, plot.ask = TRUE, ylim.cex = 1, cex.leg = 0.8, polyg = TRUE, sort.by.traits = F)
}

\arguments{
	\item{traits}{
		traits matrix with traits in column.
	}
	
	\item{var.1}{
		The first variable define the division of each plots, in most case either a vector of species or sites names.
	}
	
	\item{var.2}{
		The second variable define the division by color, in most case either a vector of species or sites names.
	}
	
	\item{col.dens}{
		A vector of colors for the second variable.
	}
	
	\item{plot.ask}{
		Logical value; ask for ploting the next plot or not.
	}
	
	\item{ylim.cex}{
		Numeric value; the magnification to be used for range of y axe
	}
	
	\item{cex.leg}{
		Numeric value; the magnification to be used for legend relative to the current setting of cex
	}
	
	\item{polyg}{
		Logical value; do the mean distribution is full or empty
	}
	
	\item{sort.by.traits}{
		Logical value; if true, sort the different plots by traits rather than by the first variable (var.1)
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
}

\examples{	
	data(finch.ind)
	plot.dens(traits.finch, sp.finch, ind.plot_finch, ylim.cex=3)
	plot.dens(traits.finch, ind.plot_finch, sp.finch, ylim.cex=8)
	plot.dens(traits.finch, rep("region", times=dim(traits.finch)[1]), sp.finch, ylim.cex=6)
	plot.dens(traits.finch, rep("toutes_sp", times=dim(traits.finch)[1]), ind.plot_finch,ylim.cex=3)
	
	plot.dens(traits.finch, sp.finch, ind.plot_finch, ylim.cex=3, sort.by.trait = TRUE)
	plot.dens(traits.finch, ind.plot_finch, sp.finch, ylim.cex=3, sort.by.trait = TRUE)
	
	plot.dens(traits.finch, ind.plot_finch, sp.finch, ylim.cex=3, sort.by.trait = FALSE)
	
	plot.dens(traits.finch, ind.plot_finch, sp.finch, ylim.cex=3, plot.ask = FALSE, polyg = FALSE)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
