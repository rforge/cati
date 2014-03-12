\name{com.index.multi}
\alias{com.index.multi}

\title{
	Calcul of custom multitraits statistics to test community assembly
}

\description{
	Calcul of custom multitraits statistics to test community assembly using null models.
}

\usage{
	com.index.multi(traits = NULL, index = NULL, by.factor = NULL, namesindex = NULL, nullmodels = NULL, ind.plot = NULL, sp = NULL, nperm = 99, printprogress = TRUE)
}

\arguments{
 
	\item{traits}{
		Individual traits matrix with traits in column.
	}
	
	\item{index}{
		A vector of function to apply to traits vectors in the form "mean(x, na.rm = TRUE)" or "range(x)".
	}
	
	\item{by.factor}{
		A factor to split the traits matrix and compute index for each subset.  
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

\seealso{
	\code{\link{com.index.c}};
	\code{\link{com.index.multi}};
	\code{\link{plot.filter}};
	\code{\link{ses}};
	\code{\link{plot.ses}}
}

\examples{

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
