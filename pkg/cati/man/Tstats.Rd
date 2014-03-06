\name{Tstats}
\alias{Tstats}
\title{
	Calcul observed T-statistics (T for Traits) and null expectations.
}
\description{
	Calcul observed T-statistics (T for Traits) as three ratios of variance, namely T_IP.IC, T_IC.IR and T_PC.PR. This function can also return the distribution of this three statistics under null models.
}

\usage{
	Tstats(Ttraits, ind_plot, sp, nperm = NULL, printprogress = TRUE, p.value = TRUE)
}

\arguments{
	\item{Ttraits}{
		Individual traits matrix with traits in column.
	}
	
	\item{ind_plot}{
		Factor defining the name of the plot in which the individu is.
	}
	
	\item{sp}{
		Factor defining the species which the individu belong to.
	}
	
	\item{nperm}{
		Number of permutation. If NULL, only observed values are returned; 
	}
	
	\item{printprogress}{
		Logical value; print progress during the calcul or not.
	}
	
	\item{p.value}{
		Logical value; calcul p-value or not.
	}
}

\details{
}

\value{
	A list of statistics:
	
	\item{$T_IP.IC}{Observed ratio between variance of individuals in populations and individuals in communities}
	\item{$T_IC.IR}{Observed ratio between variance of individuals in communities and individuals in the region}
	\item{$T_PC.PR}{Observed ratio between variance of populations in communities and individuals in the region}
	
	\item{$var_IP}{variance of individuals in populations}    
	\item{$var_PC}{variance of populations in communities}
	\item{$var_CR}{variance of communities in the region}
	\item{$var_IC}{variance of individuals in communities}
	\item{$var_PR}{variance of populations in the region}
	\item{$var_IR}{variance of individuals in the region}
	
	\item{$T_IP.IC_nm}{If nperm is numeric; Result of simulation for T_IP.IC}
	\item{$T_IC.IR_nm}{If nperm is numeric; Result of simulation for T_IC.IR}
	\item{$T_PC.PR_nm}{If nperm is numeric; Result of simulation for T_PC.PR}
	
	\item{$pval_T_IP.IC.inf}{If p.value=T; Result of p-value using the alternate hypothesis : T_IP.IC observed < T_IP.IC simulated}
	\item{$pval_T_IC.IR.inf}{If p.value=T; Result of p-value using the alternate hypothesis : T_IC.IR observed < T_IC.IR simulated}
	\item{$pval_T_PC.PR.inf}{If p.value=T; Result of p-value using the alternate hypothesis : T_PC.PR observed < T_PC.PR simulated}
	\item{$pval_T_IP.IC.sup}{If p.value=T; Result of p-value using the alternate hypothesis : T_IP.IC observed > T_IP.IC simulated}
	\item{$pval_T_IC.IR.sup}{If p.value=T; Result of p-value using the alternate hypothesis : T_IC.IR observed > T_IC.IR simulated}
	\item{$pval_T_PC.PR.sup}{If p.value=T; Result of p-value using the alternate hypothesis : T_PC.PR observed > T_PC.PR simulated}
}

\references{
	Violle, Cyrille, Brian J. Enquist, Brian J. McGill, Lin Jiang, Cecile H. Albert, Catherine Hulshof, Vincent Jung, et Julie Messier. 2012. The return of the variance: intraspecific variability in community ecology. Trends in Ecology & Evolution 27 (4): 244-252. doi:10.1016/j.tree.2011.11.014.
}

\author{
	Adrien Taudiere
	Cyrille Violle
}

\note{
}

\seealso{
	\code{\link{summary.Tstats}};
	\code{\link{plot.ses.Tstats}};
	\code{\link{plot.filter}};
	\code{\link{mat.comm}}
}

\examples{
	data(finch.ind)
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot.finch, sp=sp.finch, nperm=99)
	summary.Tstats(res.finch)
	
	#Plot functions
	plot.ses.Tstats(res.finch)
	
	# Custom theme (from rasterVis package)
	my.theme <- BuRdTheme()
	# Customize the colorkey
	my.ckey <- list(col=my.theme$regions$col)
	
	levelplot(t(ses(res.finch$T_IP.IC,res.finch$T_IP.IC_nm)$ses), colorkey=my.ckey, par.settings=my.theme,border="black")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
