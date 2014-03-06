\name{plot.filter.c}
\alias{plot.filter.c}

\title{
	Plot community assembly index for "classical" function and T-statistics
}

\description{
	Plot community assembly index and confidence intervals using T-statistics, and classical index (mean, range, CVNND & kurtosis) at the individuals and at the species levels.
}

\usage{
	plot.filter.c(com.object.ind = NULL, com.object.sp = NULL, tstats = NULL, color.cond = NULL, val.quant = c(0.025, 0.975), xlim = NULL, ylim = NULL, bysites = TRUE, cex.text = 0.8, plot.ask = TRUE, srt.text = 90)
}

\arguments{
	\item{com.object.ind}{
		The list resulting from the function com.index.c using pop.aggregate=F.
	}
	
	\item{com.object.sp}{
		The list resulting from the function com.index.c using pop.aggregate=T. 
	}
	
	\item{tstats}{
		The list resulting from the function Tstats.
	}
	
	\item{color.cond}{
		Vector of colors for conditional colors.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calcul confidence interval. By default val.quant=c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{xlim}{
		Numeric vectors of length 2, giving the x coordinates range 
	}
	
	\item{ylim}{
		Numeric vectors of length 2, giving the y coordinates range 
	}
	
	\item{bysites}{
		Logical value; plot by sites or by traits. 
	}
	
	\item{cex.text}{
		Numeric value; the magnification to be used for text relative to the current setting of cex
	}
	
	\item{plot.ask}{
		Logical value; ask for ploting the next plot or not.
	}
	
	\item{srt.text}{
		Degree of rotation for text.
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
	\code{\link{plot.filter}};
	\code{\link{plot.ses.Tstats}}
}

\examples{
data(finch.ind)
	
	#Calcul using individuals values
	res.finch.ind<-com.index.c(traits=traits.finch, ind.plot=ind.plot.finch, nperm=99, sp=sp.finch, pop.aggregate=F)
	
	#Calcul using populational values
	res.finch.sp<-com.index.c(traits=traits.finch, ind.plot=ind.plot.finch, nperm=99, sp=sp.finch, pop.aggregate=TRUE)
	
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot.finch, sp=sp.finch, nperm=99)
		
	plot.filter.c(com.object.ind=res.finch.ind, com.object.sp=res.finch.sp, tstats=res.finch, val.quant=c(0.025, 0.975), bysites=TRUE )
	plot.filter.c(com.object.ind=res.finch.ind, com.object.sp=res.finch.sp, tstats=res.finch, val.quant=c(0.025, 0.975), bysites=FALSE )
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
