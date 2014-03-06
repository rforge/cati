\name{plot.ses}
\alias{plot.ses}

\title{
	Plot standardized effect size
}

\description{
	Plot standardized effect size whith confidence interval. index.list is a list of index associate with a list of corresponding null models in this order: [1] index 1 obs - [2] index 1 null model - [3] index 2 obs - [4] index 2 null model ... Observed matrix of values need to be of the same dimension. You can transpose the observed matrix to represent either the SES by traits or by plots.
}

\usage{
	plot.ses(index.list, col.index = c("red", "purple", "green"), type = "normal", add.conf = TRUE, color.cond = TRUE, val.quant = c(0.025, 0.975), grid.v = TRUE, grid.h = TRUE, xlim = NULL, ylim = NULL)
}

\arguments{
	\item{index.list}{
		A list of index and the associate null models in the forme: list( index_1=index_1_observed, index_1_nm=null.model.index_1 ,index_2=index_2_observed, index_2_nm=null.model.index_2, ...). Each element of the list need to be named.
	}
	
	\item{col.index}{
		Vector of colors for index.
	}
	
	\item{type}{
		Type of plot. Possible type = "simple", "simple_range", "normal" and "barplot".
	}
	
	\item{add.conf}{
		Logical value; Add confidence intervals or not. 
	}
	
	\item{color.cond}{
		Logical value; If color.cond=T, color points indicate T-statistics values significatively different from the null model and grey points are not different from null model.
	}
	
	\item{val.quant}{
		Numeric vectors of length 2, giving the quantile to calcul confidence interval. By default val.quant=c(0.025,0.975) for a bilateral test with alpha = 5\%.
	}
	
	\item{grid.v}{
		Logical value; print vertical grid or not
	}
	
	\item{grid.h}{
		Logical value; print horizontal grid or not
	}
	
	\item{xlim}{
		Numeric vectors of length 2, giving the x coordinates range 
	}
	
	\item{ylim}{
		Numeric vectors of length 2, giving the y coordinates range 
	}
}

\details{
	Normal type plot T-statistics for each site, means, standard deviations, and confidence intervals of T-statistics. 
	Simple_range type plot means, standard deviations, ranges and confidence intervals of T-statistics
	Simple type type plot T-statistics for each site, means, and confidence intervals of T-statistics.
	Barplot type plot means, standard deviations and confidence intervals of T-statistics in a barplot fashion.
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
	\code{\link{plot.ses.Tstats}};
	\code{\link{plot.ses.var}}
}

\examples{
	data(finch.ind)
	res.finch<-Tstats(traits.finch, ind_plot=ind.plot.finch, sp=sp.finch, nperm=99)
	
	list.ses<-list(T_IP.IC=res.finch$T_IP.IC, T_IP.IC_nm=res.finch$T_IP.IC_nm, T_PC.PR=res.finch$T_PC.PR, T_PC.PR_nm=res.finch$T_PC.PR_nm, T_IC.IR=res.finch$T_IC.IR, T_IC.IR_nm=res.finch$T_IC.IR_nm)
	
	list.ses.t<-list(T_IP.IC=t(res.finch$T_IP.IC), T_IP.IC_nm=res.finch$T_IP.IC_nm, T_PC.PR=t(res.finch$T_PC.PR), T_PC.PR_nm=res.finch$T_PC.PR_nm, T_IC.IR=t(res.finch$T_IC.IR), T_IC.IR_nm=res.finch$T_IC.IR_nm)

	plot.ses(list.ses.t)
	
	plot.ses(list.ses)
	plot.ses(list.ses, type="simple")
	plot.ses(list.ses, type="simple_range")
	plot.ses(list.ses, type="barplot")
	
	plot.ses(list.ses, color.cond=F)
	plot.ses(list.ses, type="simple", add.conf=F)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
