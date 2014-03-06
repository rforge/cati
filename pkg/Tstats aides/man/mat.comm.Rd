\name{mat.comm}
\alias{mat.comm}

\title{
	Matrix for summarize community assembly index result
}

\description{
	Matrix for summarize community assembly index result.
}

\usage{
mat.comm(com.object = NULL, tstats = NULL, val.quant = c(0.025, 0.975), sites_names = NULL, method = "p.val")
}

\arguments{
	\item{com.object}{
		The list resulting from the function com.index.c().
	}
	
	\item{tstats}{
		The list resulting from the function Tstats.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calcul confidence interval. By default val.quant=c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{sites_names}{
		Names of sites.
	}
	
	\item{method}{
		Either "p.val" or "ses". Define the type of values in the matrix (p-value test (-1, 0 or 1) or standardized effect size).
	}
}

\details{
}

\value{
	A matrix with traits as columns and index in rows. If method="p.val", a value of 1 correspond to significatively higher value than expected under null model and a value of -1 to a significatively lower value. For an index not different from null model for a trait, the value is 0. if method="ses", the matrix is fill by standardized effect size values for each index and each traits.  
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
	\code{\link{com.index.c}};
	\code{\link{Tstats}};
	\code{\link{levelplot}}
}

\examples{
	data(finch.ind)
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot.finch, sp=sp.finch, nperm=99)
	res.finch.ind<-com.index.c(traits=traits.finch, ind.plot=ind.plot.finch, nperm=99, sp=sp.finch, pop.aggregate=F)

	mat.comm(com.object=res.finch.ind,tstats=res.finch, method="ses")
	mat.comm(com.object=res.finch.ind,tstats=res.finch, method="ses",val.quant=c(0.05,0.95))
	
	mat.comm(com.object=res.finch.ind,tstats=res.finch, method="p.val")
	mat.comm(com.object=res.finch.ind,tstats=res.finch, method="p.val",val.quant=c(0.05,0.95))
	mat.comm(com.object=res.finch.ind,tstats=res.finch, method="p.val",val.quant=c(0,0.95))
		
	require(rasterVis)
	# Custom theme (from rasterVis package)
	my.theme <- BuRdTheme()
	# Customize the colorkey
	my.ckey <- list(col=my.theme$regions$col)

	levelplot(t(mat.comm(com.object=res.finch.ind,tstats=res.finch, method="ses")), colorkey=my.ckey, par.settings=my.theme,border="black")
	levelplot(t(mat.comm(com.object=res.finch.ind,tstats=res.finch, method="p.val")), colorkey=my.ckey, par.settings=my.theme,border="black")
	levelplot(t(mat.comm(com.object=res.finch.ind,tstats=res.finch))+t(mat.comm(com.object=res.finch.ind,tstats=res.finch, val.quant=c(0.05,0.95)))+t(mat.comm(com.object=res.finch.ind,tstats=res.finch, val.quant=c(0.25,0.75))), colorkey=my.ckey, par.settings=my.theme,border="black")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
