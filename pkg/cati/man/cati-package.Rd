\name{cati-package}
\alias{cati-package}
\alias{cati}
\docType{package}

\title{
	A package for community assembly testing
}

\description{
	A package for community assembly testing in particular using individuals traits values. 
}

\details{
	\tabular{ll}{
	Package: \tab Tstats \cr
	Type: \tab Package \cr
	Version: \tab 0.1 \cr
	Date: \tab 2014-02-14 \cr
	License: \tab GPL-3 \cr
	Depends: \tab ape, nlme, e1071, moments, lattice, R(>=3.0.0) \cr
	Suggests: \tab rasterVis \cr
	}
	This package provide functions to calculate T-statistics (Tstats function) and other uni-traits metrics (com.index & com.index.c functions) to test community assembly considering the individual and other taxonomic levels. Variance partition (partvar function) and density plot for community studies (plot.dens function) are also available. Finaly, this package include function to summarise index and function to plot standardized effect size of index (e.g. plot.ses function). 
}

\author{
	Adrien Taudiere; adrien.taudiere@cefe.cnrs.fr
}

\references{
	Violle, Cyrille, Brian J. Enquist, Brian J. McGill, Lin Jiang, Cecile H. Albert, Catherine Hulshof, Vincent Jung, et Julie Messier. 2012. The return of the variance: intraspecific variability in community ecology. Trends in Ecology & Evolution 27 (4): 244-252. doi:10.1016/j.tree.2011.11.014.
}

\keyword{ package }

\seealso{
}

\examples{
	data(finch.ind)
	
	#T-statistics
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot_finch, sp=sp.finch, nperm=99)
	summary.Tstats(res.finch)
	
	list.ses<-list(T_IP.IC=res.finch$T_IP.IC, T_IP.IC_nm=res.finch$T_IP.IC_nm, T_PC.PR=res.finch$T_PC.PR, T_PC.PR_nm=res.finch$T_PC.PR_nm, T_IC.IR=res.finch$T_IC.IR, T_IC.IR_nm=res.finch$T_IC.IR_nm)
	
	plot.ses(list.ses)
	
	#Others classicals index (mean, range, CV.NND and kurtosis). Use com.index function for custom function
	res.finch.ind<-com.index.c(traits=traits.finch, ind.plot=ind.plot_finch, nperm=99, sp=sp.finch, pop.aggregate = FALSE)

	list.ses2<-list( range.ind=t(res.finch.ind$Range_obs), range.ind.null=res.finch.ind$Null_range, mean.ind=t(res.finch.ind$Mean_obs), mean.ind.null=res.finch.ind$Null_mean, CV.ind=t(res.finch.ind$CV.NND_obs), CV.ind.null=res.finch.ind$Null_CV)
	
	plot.ses(list.ses2)
	
	plot.ses(c(list.ses, list.ses2))
	plot.filter(c(list.ses, list.ses2))
}
