\name{trait.flex.anova}
\alias{trait.flex.anova}
\alias{print.trait.flex}
\alias{plot.trait.flex}

\title{
	Variance decomposition for one trait.
}

\description{
	This function decomposes the variation in community trait composition into three sources – (i) the intraspecific trait variability, (ii) the
	variability due to species turnover and (iii) their covariation is also separated. This decomposition is computed for the whole variation in the trait
	values and, The formula specified, across the contribution of various explanatory variables considered in the model. S3 method plot summarizes graphically the decomposition of trait variation, obtained with the trait.flex.anova function. Print is an other S3 method for object of class trait.flex.
}

\usage{
	trait.flex.anova(formula, specif.avg, const.avg, ...)
	\method{plot}{trait.flex}(x, plot.total= FALSE, use.percentage = TRUE, plot.covar= FALSE, cumul= FALSE, legend.pos = if (plot.total) "topleft" else "topright", plot.res = TRUE, ...)
	\method{print}{trait.flex}(x, ...)
}

\arguments{
	\item{formula}{
		The formula parameter must be a one-sided formula, i.e. starting with a tilde ( ~ ) character. The response variable is specified by the next two arguments, specif.avg and const.avg. 
	}
	
	\item{specif.avg}{
		Vector with community trait composition values for a single trait. It is calculated from trait values specific to each community (i.e. trait values for individual species are ‘specific’ to each plot, or habitat, where the species is found)
	}
	
	\item{const.avg}{
		Vector with community trait composition values for a single trait. It is calculated from average (fixed) trait values of individual species (i.e. fixed trait value for individual species used for all habitats where the species is found)
	}
	
	\item{x}{
		An object of class trait.flex.
	}
	
	\item{plot.total}{
		Logical value; if TRUE plot not only the individual components of variation, but also the total variation. This is useful particularly when the decomposition was done with non-trivial formula (i.e. with explanatory variables)
	}
	
	\item{use.percentage}{
		Logical value; if TRUE the individual plotted sources of trait variation are shown as percentages of the total variation, on 0–100 scale.
	}
	
	\item{plot.covar}{
		Logical value; if TRUE the covariance between within-species trait variability and the variability due to species composition turnover is plotted as yet another category within the stacked bars. The plot.covar argument is entirely ignored when plotting trait.flex object fitted with a formula without any predictor variables.
	}
	
	\item{cumul}{
		Logical value; if TRUE values are shown in a cumulative way.
	}
	
	\item{legend.pos}{
		This argument allows you to specify the position of graph legend. Thus argument is entirely ignored when plotting trait.flex object created with a formula without predictors
	}
	
	\item{plot.res}{
		Logical value; if resume=FALSE plot is not shown but the table of values used to print the plot is return.
	}
	
	\item{\dots}{
		Optional additional arguments.
	}
}

\details{
	The formula parameter must be a one-sided formula, i.e. starting with a tilde ( ~ ) character. The response variable is specified by the next two arguments, specif.avg and const.avg. 
}

\value{
	An object of class trait.flex. There are print and plot methods available for it. The object contains decomposition of sum of squares into intraspecific variation component, compositional variation component, their covariation and total in a SumSq element. This is a data frame with multiple rows if predictors were specified in formula argument. The RelSumSq element contains the same table relativized to unit row totals. Finally, the anova.turnover, anova.total, and anova.diff elements contain the three aov objects used to decompose the variation.
}

\references{
	Leps, Jan, Francesco de Bello, Petr Smilauer and Jiri Dolezal. 2011. Community trait response to environment: disentangling species turnover vs intraspecific trait variability effects. Ecography 34 (5): 856‑863.
}

\author{
	Jan Leps et al., 2011
}

\note{
}


\seealso{
	\code{\link{print.trait.flex}};
	\code{\link{plot.trait.flex}}
}

\examples{
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
