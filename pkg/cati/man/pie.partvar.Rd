\name{pie.partvar}
\alias{pie.partvar}

\title{
	Pie plot of variance partitionning.
}

\description{
	Plot pies of the partition of variance calculate with the function partvar.
}

\usage{
	pie.partvar(var.part, col.pie = NA)
}

\arguments{
	\item{var.part}{
	The result of the partvar function
	}
	\item{col.pie}{
	Vector of color for pie
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

\seealso{
	\code{\link{partvar}};
	\code{\link{bar.partvar}}
}

\examples{
	data(finch.ind)
	res.partvar.finch<-partvar(traits=traits.finch, factors=cbind(as.factor(as.vector(plot.partvar)), as.factor(as.vector(sp.partvar))))
	pie.partvar(res.partvar.finch, col=c("red", "green", "blue"))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
