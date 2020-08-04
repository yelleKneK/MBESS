\name{ss.aipe.crd}
\alias{ss.aipe.crd.nclus.fixedwidth}
\alias{ss.aipe.crd.nindiv.fixedwidth}
\alias{ss.aipe.crd.nclus.fixedbudget}
\alias{ss.aipe.crd.nindiv.fixedbudget}
\alias{ss.aipe.crd.both.fixedbudget}
\alias{ss.aipe.crd.both.fixedwidth}
\title{
	Find target sample sizes for the accuracy in unstandardized conditions means estimation in CRD
}
\description{
	Find target sample sizes (the number of clusters, cluster size, or both) for the accuracy in unstandardized conditions means estimation in CRD. If users wish to seek for both types of sample sizes simultaneously, an additional constraint is required, such as a desired width or a desired budget.
}
\usage{
ss.aipe.crd.nclus.fixedwidth(width, nindiv, prtreat, tauy=NULL, sigma2y=NULL, 
	totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, 
	assurance=NULL, conf.level = 0.95, cluscost=NULL, indivcost=NULL, diffsize=NULL)
ss.aipe.crd.nindiv.fixedwidth(width, nclus, prtreat, tauy=NULL, sigma2y=NULL, 
	totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, 
	assurance=NULL, conf.level = 0.95, cluscost=NULL, indivcost=NULL, diffsize=NULL)
ss.aipe.crd.nclus.fixedbudget(budget, nindiv, cluscost = 0, indivcost = 1, 
	prtreat = NULL, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, 
	r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, diffsize=NULL)
ss.aipe.crd.nindiv.fixedbudget(budget, nclus, cluscost = 0, indivcost = 1, 
	prtreat = NULL, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, 
	r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, diffsize=NULL)
ss.aipe.crd.both.fixedbudget(budget, cluscost=0, indivcost=1, prtreat, tauy=NULL, 
	sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, 
	numpredictor = 0, assurance=NULL, conf.level = 0.95, diffsize=NULL)
ss.aipe.crd.both.fixedwidth(width, cluscost=0, indivcost=1, prtreat, tauy=NULL, 
	sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, 
	numpredictor = 0, assurance=NULL, conf.level = 0.95, diffsize=NULL)
}
\arguments{
  \item{width}{
	The desired width of the confidence interval of the unstandardized means difference
  }
  \item{budget}{
	The desired amount of budget
  }
  \item{nclus}{
	The desired number of clusters
  }
  \item{nindiv}{
	The number of individuals in each cluster (cluster size)
  }
  \item{prtreat}{
	The proportion of treatment clusters
  }
  \item{cluscost}{
	The cost of collecting a new cluster regardless of the number of individuals collected in each cluster
  }
  \item{indivcost}{
	The cost of collecting a new individual
  }
  \item{tauy}{
	The residual variance in the between level before accounting for the covariate
  }
  \item{sigma2y}{
	The residual variance in the within level before accounting for the covariate
  }
  \item{totalvar}{
	The total resiudal variance before accounting for the covariate
  }
  \item{iccy}{
	The intraclass correlation of the dependent variable
  }
  \item{r2within}{
	The proportion of variance explained in the within level (used when \code{covariate = TRUE})
  }
  \item{r2between}{
	The proportion of variance explained in the between level (used when \code{covariate = TRUE})
  }
  \item{numpredictor}{
	The number of predictors used in the between level
  }
  \item{assurance}{
	The degree of assurance, which is the value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g, .80, .90, .95)
  }
  \item{conf.level}{
	The desired level of confidence for the confidence interval
  }
  \item{diffsize}{
	Difference cluster size specification. The difference in cluster sizes can be specified in two ways. First, users may specify cluster size as integers, which can be negative or positive. The resulting cluster sizes will be based on the estimated cluster size adding by the specified vectors. For example, if the cluster size is 25, the number of clusters is 10, and the specified different cluster size is \code{c(-1, 0, 1)}, the cluster sizes will be 24, 25, 26, 24, 25, 26, 24, 25, 26, and 24. Second, users may specify cluster size as positive decimals. The resulting cluster size will be based on the estimated cluster size multiplied by the specified vectors. For example, if the cluster size is 25, the number of clusters is 10, and the specified different cluster size is \code{c(-1, 0, 1)}, the cluster sizes will be 24, 25, 26, 24, 25, 26, 24, 25, 26, and 24. If \code{NULL}, the cluster size is equal across clusters. 
  }
}
\details{
	Here are the functions' descriptions:
	\itemize{
		\item \code{ss.aipe.crd.nclus.fixedwidth} Find the number of clusters given a specified width of the confidence interval and the cluster size
		\item \code{ss.aipe.crd.nindiv.fixedwidth} Find the cluster size given a specified width of the confidence interval and the number of clusters
		\item \code{ss.aipe.crd.nclus.fixedbudget} Find the number of clusters given a budget and the cluster size
		\item \code{ss.aipe.crd.nindiv.fixedbudget} Find the cluster size given a budget and the number of clusters
		\item \code{ss.aipe.crd.both.fixedbudget} Find the sample size combinations (the number of clusters and that cluster size) providing the narrowest confidence interval given the fixed budget
		\item \code{ss.aipe.crd.both.fixedwidth} Find the sample size combinations (the number of clusters and that cluster size) providing the lowest cost given the specified width of the confidence interval
	}
}
\value{
	The \code{ss.aipe.crd.nclus.fixedwidth} and \code{ss.aipe.crd.nclus.fixedbudget} functions provide the number of clusters. The \code{ss.aipe.crd.nindiv.fixedwidth} and \code{ss.aipe.crd.nindiv.fixedbudget} functions provide the cluster size. The \code{ss.aipe.crd.both.fixedbudget} and \code{ss.aipe.crd.both.fixedwidth} provide the number of clusters and the cluster size, respectively.
}
\references{
Pornprasertmanic, S., & Schneider, W. J. (2014). Accuracy in parameter estimation in cluster randomized designs. \emph{Psychological Methods}, \emph{19}, 356--379.

}
\author{
Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
\dontrun{
# Examples for each function
ss.aipe.crd.nclus.fixedwidth(width=0.3, nindiv=30, prtreat=0.5, tauy=0.25, sigma2y=0.75)

ss.aipe.crd.nindiv.fixedwidth(width=0.3, nclus=250, prtreat=0.5, tauy=0.25, sigma2y=0.75)

ss.aipe.crd.nclus.fixedbudget(budget=10000, nindiv=20, cluscost=20, indivcost=1)

ss.aipe.crd.nindiv.fixedbudget(budget=10000, nclus=30, cluscost=20, indivcost=1, 
prtreat=0.5, tauy=0.05, sigma2y=0.95, assurance=0.8)

ss.aipe.crd.both.fixedbudget(budget=10000, cluscost=30, indivcost=1, prtreat=0.5, tauy=0.25, 
	sigma2y=0.75)

ss.aipe.crd.both.fixedwidth(width=0.3, cluscost=0, indivcost=1, prtreat=0.5, tauy=0.25, 
	sigma2y=0.75)

# Examples for different cluster size
ss.aipe.crd.nclus.fixedwidth(width=0.3, nindiv=30, prtreat=0.5, tauy=0.25, sigma2y=0.75, 
diffsize = c(-2, 1, 0, 2, -1, 3, -3, 0, 0))

ss.aipe.crd.nclus.fixedwidth(width=0.3, nindiv=30, prtreat=0.5, tauy=0.25, sigma2y=0.75, 
diffsize = c(0.6, 1.2, 0.8, 1.4, 1, 1, 1.1, 0.9))
}
}