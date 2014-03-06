\name{com.index}
\alias{com.index}

\title{
	Calcul of custom statistics to test community assembly
}

\description{
	Calcul of custom statistics to test community assembly using null models.
}

\usage{
	com.index(traits = NULL, index = NULL, namesindex = NULL, nullmodels = NULL, ind.plot = NULL, sp = NULL, nperm = 99, printprogress = TRUE)
}

\arguments{
	
	\item{traits}{
		Individual traits matrix with traits in column.
	}
	
	\item{index}{
		A vector of function to apply to traits vectors in the form "mean(x, na.rm=T)" or "range(x)".
	}
	
	\item{namesindex}{
		A vector of names for index.
	}
	
	\item{nullmodels}{
		A vector of values for null models coresponding to index. A value of 1 correspond to randomization of individual values whithin community. A value of 2 correspond to randomization of individual values whithin region. A values of 3 correspond to randomization of population values whithin region.
	}
	
	\item{ind.plot}{
		Factor defining the name of the plot in which the individu is. One of comm and ind_plot is sufficient.
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
}

\details{
Calcul of statistics (e.g. mean, range, CVNND and kurtosis) to test community assembly using null models.
For each statistic this function return observed value and correspondant Null distribution.
This function implement three null models which keep unchanged the number of individual per community.
Models 1 correspond to randomization of individual values whithin community.
Models 2 correspond to randomization of individual values whithin region.
Models 3 correspond to randomization of population values whithin region.

In most case, model 1 and 2 correspond to index at the individual level and the model 3 to index at the species (or any other aggregate variable like genus or family) level.
}

\value{
	A list of list:
	
	\item{$obs}{List of observed values for each traits in each community. Each component of the list correspond to a matrix containing the result for each custom function.}
	
	\item{$Null}{List of null values for each traits in each community. Each component of the list correspond to an array containing the result of the permutation for each custom function.}
	
	\item{$list.index}{List of index values and associate null models. Internal use in other function. Traits in columns.}
	\item{$list.index.t}{List of index values and associate null models. Internal use in other function. Traits in rows.}
	
	\item{$sites_richness}{Number of species by sites.}
	\item{$namestraits}{Names of traits.}
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
	\code{\link{plot.filter}};
	\code{\link{ses}};
	\code{\link{plot.ses}}
}

\examples{
	data(finch.ind)
	
	#_________________
	#simple function (mean and var)
	#Null model 1 is trivial for this simple function because randomisation is whithin community only
	
	res.finch.ind_mean_var1<-com.index(traits=traits.finch, index=c("mean(x, na.rm=T)", "var(x, na.rm=T)"), sp=sp.finch , nullmodels=c(1,1), ind.plot=ind.plot.finch, nperm=NULL)
	res.finch.ind_mean_var1$obs$mean
	
	res.finch.ind_mean_var2<-com.index(traits=traits.finch, index=c("mean(x, na.rm=T)", "var(x, na.rm=T)"), sp=sp.finch , nullmodels=c(2,2), ind.plot=ind.plot.finch, nperm=99)
	res.finch.ind_mean_var2$obs$mean
	res.finch.ind_mean_var2$Null$mean
	ses(res.finch.ind_mean_var2$obs$mean, res.finch.ind_mean_var2$Null$mean)
	
	
	res.finch.ind_mean_var3<-com.index(traits=traits.finch, index=c("mean(x, na.rm=T)", "var(x, na.rm=T)"), sp=sp.finch , nullmodels=c(3,3), ind.plot=ind.plot.finch, nperm=99)
	res.finch.ind_mean_var3$obs$mean
	ses(res.finch.ind_mean_var3$obs$mean, res.finch.ind_mean_var3$Null$mean)
	
	#_________________
	#More complex function: mean & var by population
	
	res.finch.ind_mean_var.by.pop1<-com.index(traits=traits.finch, index=c("tapply(x, paste(sp.finch, ind.plot.finch), function(x) mean(x, na.rm=T))", "tapply(x, paste(sp.finch, ind.plot.finch), function(x) sd(x, na.rm=T))"), sp=sp.finch , nullmodels=c(1,1), ind.plot=ind.plot.finch, nperm=99, namesindex=c("mean_1","sd_1"))
	ses(res.finch.ind_mean_var.by.pop1$obs$mean, res.finch.ind_mean_var.by.pop1$Null$mean)$ses
	plot.ses(res.finch.ind_mean_var.by.pop1$list.index)
	
	res.finch.ind_mean_var.by.pop2<-com.index(traits=traits.finch, index=c("tapply(x, paste(sp.finch, ind.plot.finch), function(x) mean(x, na.rm=T))", "tapply(x, paste(sp.finch, ind.plot.finch), function(x) sd(x, na.rm=T))"), sp=sp.finch , nullmodels=c(2,2), ind.plot=ind.plot.finch, nperm=99, namesindex=c("mean_2","sd_2"))
	ses(res.finch.ind_mean_var.by.pop2$obs$mean, res.finch.ind_mean_var.by.pop2$Null$mean)$ses
	plot.ses(res.finch.ind_mean_var.by.pop2$list.index)
	
	plot.ses(c(res.finch.ind_mean_var.by.pop1$list.index, res.finch.ind_mean_var.by.pop2$list.index))
	
	#_________________
	#More complex function: mean by sites using population means
	
	sites.pop<-unlist(strsplit(rownames(res.finch.ind_mean_var.by.pop1$obs$mean), split=" "))[seq(2,2*length(rownames(res.finch.ind_mean_var.by.pop1$obs$mean)),by=2)]
	apply(res.finch.ind_mean_var.by.pop1$obs$mean, 2, function(x) tapply(x, sites.pop, mean, na.rm=T)  )
	
	ses(apply(res.finch.ind_mean_var.by.pop1$obs$mean, 2, function(x) tapply(x, sites.pop, mean, na.rm=T)), apply(res.finch.ind_mean_var.by.pop1$Null$mean, c(1,3), function(x) tapply(x, sites.pop, mean, na.rm=T)  ))
	
	plot.ses(list(Mean.by.pop_obs=apply(res.finch.ind_mean_var.by.pop1$obs$mean, 2, function(x) tapply(x, sites.pop, mean, na.rm=T)), Mean.by.pop_Null=apply(res.finch.ind_mean_var.by.pop1$Null$mean, c(1,3), function(x) tapply(x, sites.pop, mean, na.rm=T)  ), sd.by.pop_obs=apply(res.finch.ind_mean_var.by.pop1$obs$sd, 2, function(x) tapply(x, sites.pop, mean, na.rm=T)), sd.by.pop_Null=apply(res.finch.ind_mean_var.by.pop1$Null$sd, c(1,3), function(x) tapply(x, sites.pop, mean, na.rm=T)  )  ) )

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
