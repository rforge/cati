\name{plot.filter}
\alias{plot.filter}

\title{
	Plot community assembly index for custom functions
}

\description{
	Plot community assembly index and confidence intervals using a list of index.
}

\usage{
	plot.filter(index.list, color.cond = NULL, val.quant = c(0.025, 0.975), cex.text = 0.8, plot.ask = TRUE, srt.text = 90, xlim = NULL, ylim = NULL, bysites = T)
}

\arguments{
	
	\item{index.list}{
		A list of index and the associate null models in the forme: list( index_1=index_1_observed, index_1_nm=null.model.index_1 ,index_2=index_2_observed, index_2_nm=null.model.index_2, ...). 
	}
	
	\item{color.cond}{
		Vector of colors for conditional colors.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calcul confidence interval. By default val.quant=c(0.025,0.975) for a bilateral test with alpha = 5\%.
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
	
	\item{xlim}{
		Numeric vectors of length 2, giving the x coordinates range
	}
	
	\item{ylim}{
		Numeric vectors of length 2, giving the y coordinates range
	}
	
	\item{bysites}{
		Logical value; plot by sites or by traits. 
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
	\code{\link{plot.ses.Tstats}}
}

\examples{	
	data(finch.ind)
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot_finch, sp=sp.finch, nperm=99)
	res.finch.sp<-com.index.c(traits = traits.finch, ind.plot = ind.plot_finch, nperm = 99, sp = sp.finch, pop.aggregate = TRUE)
	res.finch.ind<-com.index.c(traits = traits.finch, ind.plot = ind.plot_finch, nperm = 99, sp = sp.finch, pop.aggregate = FALSE)
	
	index.list<-list(T_IP.IC=res.finch$T_IP.IC, T_IP.IC_nm=res.finch$T_IP.IC_nm, 
	                 T_IC.IR=res.finch$T_IC.IR, T_IC.IR_nm=res.finch$T_IC.IR_nm, 
	                 T_PC.PR=res.finch$T_PC.PR, T_PC.PR_nm=res.finch$T_PC.PR_nm,
	                 range.ind=t(res.finch.ind$Range_obs), range.ind.null=res.finch.ind$Null_range,
	                 mean.ind=t(res.finch.ind$Mean_obs), mean.ind.null=res.finch.ind$Null_mean,
	                 CV.ind=t(res.finch.ind$CV.NND_obs), CV.ind.null=res.finch.ind$Null_CV,
	                 kurtosis.ind=t(res.finch.ind$kurtosis_obs), kurtosis.ind.null=res.finch.ind$Null_kurtosis,
	                 range.sp=t(res.finch.sp$Range_obs), range.sp.null=res.finch.sp$Null_range,
	                 mean.sp=t(res.finch.sp$Mean_obs), mean.sp.null=res.finch.sp$Null_mean,
	                 CV.sp=t(res.finch.sp$CV.NND_obs), CV.sp.null=res.finch.sp$Null_CV,
	                 kurtosis.sp=t(res.finch.sp$kurtosis_obs), kurtosis.sp.null=res.finch.sp$Null_kurtosis)         
	
	plot.filter(index.list)
	
	plot.filter(list(T_IP.IC=res.finch$T_IP.IC, T_IP.IC_nm=res.finch$T_IP.IC_nm, 
	                 T_IC.IR=res.finch$T_IC.IR, T_IC.IR_nm=res.finch$T_IC.IR_nm, 
	                 T_PC.PR=res.finch$T_PC.PR, T_PC.PR_nm=res.finch$T_PC.PR_nm))
	
	plot.filter(list(T_IP.IC=res.finch$T_IP.IC, T_IP.IC_nm=res.finch$T_IP.IC_nm, 
	                 T_IC.IR=res.finch$T_IC.IR, T_IC.IR_nm=res.finch$T_IC.IR_nm, 
	                 T_PC.PR=res.finch$T_PC.PR, T_PC.PR_nm=res.finch$T_PC.PR_nm), bysites=FALSE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
