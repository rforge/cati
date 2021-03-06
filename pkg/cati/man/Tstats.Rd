\name{Tstats}
\alias{Tstats}
\alias{barplot.Tstats}
\alias{plot.Tstats}
\alias{summary.Tstats}

\title{
	Calcul observed T-statistics (T for Traits) and null expectations.
}
\description{
	Calcul observed T-statistics (T for Traits) as three ratios of variance, namely T_IP.IC, T_IC.IR and T_PC.PR. This function can also return the distribution of this three statistics under null models.
}

\usage{
	Tstats(Ttraits, ind_plot, sp, nperm = NULL, printprogress = TRUE, p.value = TRUE)
	
	\method{barplot}{Tstats}(height, val.quant = c(0.025, 0.975), col.Tstats = c("red", "purple", "green", "white"), ylim=NULL, \dots)
	\method{plot}{Tstats}(x, val.quant = c(0.025, 0.975), col.Tstats = c("red", "purple", "green"), type = "normal", add.conf = TRUE, \dots)
	\method{summary}{Tstats}(object, val.quant = c(0.025, 0.975), type = "all")

}

\arguments{
	\item{Ttraits}{
		Individual Traits matrix with traits in columns.
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
	
	\item{x}{
		An object of class Tstats.
	}
	
	\item{height}{
		An object of class Tstats.
	}
	
	\item{object}{
		An object of class Tstats.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calcul confidence interval. By default val.quant=c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{col.Tstats}{
		A vector of three color correspond to the three T-statistics.
	}
	
	\item{type}{
		For the plot function, type of plot. Either "color_cond", "simple", "simple_sd", "normal" and "barplot". 
		For the summary function, type of summary statistics. Either "binary", "percent", "p.value", "site" or "all". 
	}
	
	\item{add.conf}{
		Logical value; Add confidence intervals or not. 
	}
	
	\item{ylim}{
		limits for the y axis.
	}
	
	\item{\dots}{
		Any additional arguments are passed to the plot function creating the core of the plot and can be used to adjust the look of resulting graph.
	}
}

\details{
	S3 method plot:
	-Normal type plot means, standard deviations, ranges and confidence intervals of T-statistics.
	-Color_cond type plot T-statistics for each site and traits with color for significative values and grey for non signifivative ones. 
	-Simple_sd type plot means, standard deviations and confidence intervals of T-statistics
	-Simple type plot T-statistics for each site and traits and the mean confidence intervals by traits
	-Barplot type plot means, standard deviations and confidence intervals of T-statistics in a barplot fashion
	
	S3 method summary:
	-Binary type only test if a T-statistics is significatively different from the null expectation for each traits.
	-Percent type determine the percentage of sites were the T-statistics is significatively different from the null expectation for each traits. Asterix shows global significance of the test.
	-P-value type determine the p-value (two unilateral tests) of the T-statistics for each traits and sites.
	-Site type allows to know in which sites T-statistics deviate from the null expectation.
	-All type do all the precedent type of summary.
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
	
	\item{$pval_T_IP.IC.inf}{If p.value=TRUE; Result of p-value using the alternate hypothesis : T_IP.IC observed < T_IP.IC simulated}
	\item{$pval_T_IC.IR.inf}{If p.value=TRUE; Result of p-value using the alternate hypothesis : T_IC.IR observed < T_IC.IR simulated}
	\item{$pval_T_PC.PR.inf}{If p.value=TRUE; Result of p-value using the alternate hypothesis : T_PC.PR observed < T_PC.PR simulated}
	\item{$pval_T_IP.IC.sup}{If p.value=TRUE; Result of p-value using the alternate hypothesis : T_IP.IC observed > T_IP.IC simulated}
	\item{$pval_T_IC.IR.sup}{If p.value=TRUE; Result of p-value using the alternate hypothesis : T_IC.IR observed > T_IC.IR simulated}
	\item{$pval_T_PC.PR.sup}{If p.value=TRUE; Result of p-value using the alternate hypothesis : T_PC.PR observed > T_PC.PR simulated}
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
	\code{\link{com.index}};
	\code{\link{com.index.multi}};
	\code{\link{plot_cor.Tstats}};
	\code{\link{plot_ses.var}};
	\code{\link{plot.listofindex}}
}

\examples{
	data(finch.ind)
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot.finch, sp=sp.finch, nperm=9, print=FALSE)
	attributes(res.finch)

	#Tstats class is associated to S3 methods plot, barplot and summary
	
	plot(res.finch)
	
	plot(res.finch, type="color_cond")
	plot(res.finch, type="simple")
	plot(res.finch, type="simple_sd")
	plot(res.finch, type="barplot")
	
	attributes(summary(res.finch))
	head(summary(res.finch)$p.value, 10)
	
	summary.Tstats(res.finch, type="binary")
	summary.Tstats(res.finch, type="percent")
	summary.Tstats(res.finch, type="site")
	summary.Tstats(res.finch, type="p.value")
	summary.Tstats(res.finch, type="all")
	
	barplot(res.finch)
	
	attributes(summary(res.finch))
	head(summary(res.finch)$p.value, 10)
	
	#### An other way to see "ses values" of T-statistics
	
	# Custom theme (from rasterVis package)
	my.theme <- BuRdTheme()
	# Customize the colorkey
	my.ckey <- list(col=my.theme$regions$col)
	
	levelplot(t(ses(res.finch$T_IP.IC,res.finch$T_IP.IC_nm)$ses), colorkey=my.ckey, par.settings=my.theme,border="black")
}



\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
