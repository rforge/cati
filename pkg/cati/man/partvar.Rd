\name{partvar}
\alias{partvar}

\title{
	Variance partitioning accross nested scales
}

\description{
	Variance partitioning accross nested scales using a decomposition (varcomp function) of variance on restricted maximum likelihood (REML) method (lme function). See Messier et al. 2010 for more information. 
}

\usage{
	partvar(traits, factors, printprogress = TRUE)
}

\arguments{
	\item{traits}{
		traits matrix with traits in column
	}
  
	\item{factors}{
		A matrix of factors with the first column correspond to the higher level factor, the second row the second higher level factor and so on.  
	}
	
	\item{printprogress}{
		Logical value; print progress during the calcul or not.
	}
}

\details{
}

\value{
	A matrix of variance values with traits in rows and nested factors in column 
}

\references{
	Messier, Julie, Brian J. McGill, et Martin J. Lechowicz. 2010. How do traits vary across ecological scales? A case for trait-based ecology: How do traits vary across ecological scales? Ecology Letters 13(7): 838-848. doi:10.1111/j.1461-0248.2010.01476.x.
}

\author{
	Adrien Taudiere
	Julie Messier
}

\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{pie.partvar}};
	\code{\link{bar.partvar}}
}

\examples{
	data(finch.ind)
	res.partvar.finch<-partvar( traits=traits.finch, factors=cbind(as.factor(as.vector(ind.plot_finch)), as.factor(as.vector(sp.finch)), as.factor(as.vector(unlist(strsplit(as.vector(sp.finch),"_")) [seq(1,length(sp.finch)*2, by=2)])) ) )
	
	pie.partvar(res.partvar.finch, col=c("red", "green", "blue"))

	bar.partvar(res.partvar.finch, col=c("red", "green", "blue"))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
