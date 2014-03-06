\name{plot.com}
\alias{plot.com}

\title{
	Plot community assembly index for the "classical" function
}

\description{
	Plot community assembly index for the "classical" function com.index.c (mean, range, CVNND and kurtosis). If names of sites are NULL, the label one the plot correspond to the site's species richness.
}

\usage{
	plot.com(com.object = NULL, namestraits = NULL, sites_names = NULL, plot.ask = TRUE, ord = NULL)
}

\arguments{
	\item{com.object}{
		The list resulting from the function com.index.c using pop.aggregate=F.
	}
	
	\item{namestraits}{
		Names of traits.
	}
	
	\item{sites_names}{
		Names of traits.
	}
	
	\item{plot.ask}{		
		Logical value; ask for ploting the next plot or not.
	}
	
	\item{ord}{
		A vector to order sites.
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
	\code{\link{plot.filter.c}};
	\code{\link{plot.ses}}
}

\examples{
	data(finch.ind)
	res.finch.ind<-com.index.c(traits=traits.finch, ind.plot=ind.plot.finch, nperm=99, sp=sp.finch, pop.aggregate=F)
	plot.com(res.finch.ind)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
