\name{com.index.c}
\alias{com.index.c}

\title{
	Classical index of community assembly
}

\description{
	Calcul of four statistics (mean, range, CVNND and kurtosis) to test community assembly using null models. For each statistic this function calcul the mean, the standard deviation and two quantiles.
}

\usage{
	com.index.c(traits = NULL, ind.plot = NULL, sp = NULL, proba = NULL, nperm = 99, val.quant = c(0.025, 0.975), pop.aggregate = TRUE, printprogress = TRUE)
}

\arguments{
	\item{traits}{
		Individual traits matrix with traits in column.
	}
	
	\item{sp}{
		Factor defining the species which the individu belong to.
	}
	
	\item{ind.plot}{
		Factor defining the name of the plot in which the individu is. One of comm and ind_plot is sufficient.
	}
	
	\item{proba}{
		Proba can be use to take into account for abundance (c.f. Jung et al. 2010)
	}
	
	\item{nperm}{
		Number of permutation. If NULL, only observed values are returned; 
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calcul confidence interval. By default val.quant=c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{pop.aggregate}{
		Logical value; Calcul the index using population traits (mean of all individuals of one species in one site) instead of individuals traits. 
	}
	
	\item{printprogress}{
		Logical value; print progress during the calcul or not.
	}
}

\details{
}

\value{
	A list of statistics:
	
	\item{$Range_obs}{Observed range values for each traits in each community.}
	\item{$Mean_obs}{Observed mean values for each traits in each community.}
	\item{$CV.NND_obs}{Observed CVNND values for each traits in each community.}
	\item{$kurtosis_obs}{Observed kurtosis values for each traits in each community.}
	
	\item{$Range_Mean}{Mean of null range values for each traits in each community.}
	\item{$Mean_Mean}{Mean of null mean values for each traits in each community.}
	\item{$CV.NND_Mean}{Mean of null CV.NND values for each traits in each community.}
	\item{$kurtosis_Mean}{Mean of null kurtosis values for each traits in each community.}

	\item{$Range_sd}{Standard deviation of null range values for each traits in each community.}
	\item{$Mean_sd}{Standard deviation of null mean values for each traits in each community.}
	\item{$CV.NND_sd}{Standard deviation of null CV.NND values for each traits in each community.}
	\item{$kurtosis_sd}{Standard deviation of null kurtosis values for each traits in each community.}
	
	\item{$Range_qsup}{Upper quantile of null range values for each traits in each community.}
	\item{$Mean_qsup}{Upper quantile of null mean values for each traits in each community.}
	\item{$CV.NND_qsup}{Upper quantile of null CV.NND values for each traits in each community.}
	\item{$kurtosis_qsup}{Upper quantile of null kurtosis values for each traits in each community.}
	
	\item{$Range_qinf}{Lower quantile of null range values for each traits in each community.}
	\item{$Mean_qinf}{Lower quantile of null mean values for each traits in each community.}
	\item{$CV.NND_qinf}{Lower quantile of null CV.NND values for each traits in each community.}
	\item{$kurtosis_qinf}{Lower quantile of null kurtosis values for each traits in each community.}
	
	\item{$Null_range}{A list of permutation of range values for each traits in each community.}
	\item{$Null_mean}{A list of permutation of mean values for each traits in each community.}
	\item{$Null_CV}{A list of permutation of CVNND values for each traits in each community.}
	\item{$Null_kurtosis}{A list of permutation of kurtosis values for each traits in each community.}
	
	\item{$sites_richness}{Number of species by sites}.
	\item{$namestraits}{Names of traits.}
}

\references{
	Jung, Vincent, Cyrille Violle, Cedric Mondy, Lucien Hoffmann, et Serge Muller. 2010. Intraspecific variability and trait-based community assembly: Intraspecific variability and community assembly. Journal of Ecology 98 (5): 1134-1140. doi:10.1111/j.1365-2745.2010.01687.x.
}

\author{
	Adrien Taudiere
	Vincent Jung
}

\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{com.index}};
	\code{\link{plot.com}};
	\code{\link{plot.filter.c}}
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
	
	plot.com(res.finch.ind)	
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
