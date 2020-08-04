\name{ss.aipe.crd.es}
\alias{ss.aipe.crd.es.nclus.fixedwidth}
\alias{ss.aipe.crd.es.nindiv.fixedwidth}
\alias{ss.aipe.crd.es.nclus.fixedbudget}
\alias{ss.aipe.crd.es.nindiv.fixedbudget}
\alias{ss.aipe.crd.es.both.fixedbudget}
\alias{ss.aipe.crd.es.both.fixedwidth}
\title{
	Find target sample sizes for the accuracy in standardized conditions means estimation in CRD
}
\description{
	Find target sample sizes (the number of clusters, cluster size, or both) for the accuracy in standardized conditions means estimation in CRD. If users wish to seek for both types of sample sizes simultaneously, an additional constraint is required, such as a desired width or a desired budget. This function uses the likelihood-based confidence interval (Cheung, 2009) by the \code{OpenMx} package (Boker et al., 2011). See further details at Pornprasertmanit and Schneider (2010, submitted).
}
\usage{
ss.aipe.crd.es.nclus.fixedwidth(width, nindiv, es, estype = 1, iccy, prtreat, 
	r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, 
	conf.level = 0.95, nrep = 1000, iccz = NULL, seed = 123321, multicore = FALSE, 
	numProc=NULL, cluscost=NULL, indivcost=NULL, diffsize=NULL)
ss.aipe.crd.es.nindiv.fixedwidth(width, nclus, es, estype = 1, iccy, prtreat, 
	r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, 
	conf.level = 0.95, nrep = 1000, iccz = NULL, seed = 123321, multicore = FALSE, 
	numProc=NULL, cluscost=NULL, indivcost=NULL, diffsize=NULL) 
ss.aipe.crd.es.nclus.fixedbudget(budget, nindiv, cluscost, indivcost, nrep=NULL, 
	prtreat=NULL, iccy=NULL, es=NULL, estype = 1, numpredictor = 0, 
	iccz=NULL, r2within=NULL, r2between=NULL, assurance=NULL, 
	seed=123321, multicore=FALSE, numProc=NULL, conf.level=0.95, diffsize=NULL)
ss.aipe.crd.es.nindiv.fixedbudget(budget, nclus, cluscost, indivcost, nrep=NULL, 
	prtreat=NULL, iccy=NULL, es=NULL, estype = 1, numpredictor = 0, 
	iccz=NULL, r2within=NULL, r2between=NULL, assurance=NULL, 
	seed=123321, multicore=FALSE, numProc=NULL, conf.level=0.95, diffsize=NULL) 
ss.aipe.crd.es.both.fixedbudget(budget, cluscost=0, indivcost=1, es, estype = 1, 
	iccy, prtreat, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, 
	conf.level = 0.95, nrep = 1000, iccz = NULL, seed = 123321, multicore = FALSE, 
	numProc=NULL, diffsize=NULL)
ss.aipe.crd.es.both.fixedwidth(width, cluscost=0, indivcost=1, es, estype = 1, 
	iccy, prtreat, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, 
	conf.level = 0.95, nrep = 1000, iccz = NULL, seed = 123321, multicore = FALSE, 
	numProc=NULL, diffsize=NULL)
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
  \item{iccy}{
	The intraclass correlation of the dependent variable
  }
  \item{es}{
	The amount of effect size
  }
  \item{estype}{
	The type of effect size. There are only three possible options: 0 = the effect size using total standard deviation, 1 = the effect size using the individual-level standard deviation (level 1), 2 = the effect size using the cluster-level standard deviation (level 2)
  }
  \item{numpredictor}{
	If 1, a single covariate is included into the model. If 0, the no-covariate model is used. This function cannot handle multiple covariates. Therefore, only the values of 0 and 1 are allowed.
  }
  \item{iccz}{
	The intraclass correlation of the covariate (used when \code{covariate = TRUE}). If \code{iccz = 0}, the within-level covariate will be only used. If \code{iccz = 1}, the between-level covariate will be only used.
  }
  \item{r2within}{
	The proportion of variance explained in the within level (used when \code{covariate = TRUE})
  }
  \item{r2between}{
	The proportion of variance explained in the between level (used when \code{covariate = TRUE})
  }
  \item{assurance}{
	The degree of assurance, which is the value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g, .80, .90, .95)
  }
  \item{nrep}{
	The number of replications used in a priori Monte Carlo simulation
  }
  \item{seed}{
	A desired seed number
  }
  \item{multicore}{
	Use multiple processors within a computer. Specify as \code{TRUE} to use it. 
  }
  \item{numProc}{
	The number of processors to be used when \code{multicore=TRUE}. If it is not specified, the package will use the maximum number of processors in a machine. 
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
		\item \code{ss.aipe.crd.es.nclus.fixedwidth} Find the number of clusters given a specified width of the confidence interval and the cluster size
		\item \code{ss.aipe.crd.es.nindiv.fixedwidth} Find the cluster size given a specified width of the confidence interval and the number of clusters
		\item \code{ss.aipe.crd.es.nclus.fixedbudget} Find the number of clusters given a budget and the cluster size
		\item \code{ss.aipe.crd.es.nindiv.fixedbudget} Find the cluster size given a budget and the number of clusters
		\item \code{ss.aipe.crd.es.both.fixedbudget} Find the sample size combinations (the number of clusters and that cluster size) providing the narrowest confidence interval given the fixed budget
		\item \code{ss.aipe.crd.es.both.fixedwidth} Find the sample size combinations (the number of clusters and that cluster size) providing the lowest cost given the specified width of the confidence interval
	}
}
\value{
	The \code{ss.aipe.crd.es.nclus.fixedwidth} and \code{ss.aipe.crd.es.nclus.fixedbudget} functions provide the number of clusters. The \code{ss.aipe.crd.es.nindiv.fixedwidth} and \code{ss.aipe.crd.es.nindiv.fixedbudget} functions provide the cluster size. The \code{ss.aipe.crd.es.both.fixedbudget} and \code{ss.aipe.crd.es.both.fixedwidth} provide the number of clusters and the cluster size, respectively.
}
\references{ 
Boker, S., M., N., Maes, H., Wilde, M., Spiegel, M., Brick, T., et al. (2011).
OpenMx: An open source extended structural equation modeling framework.
\emph{Psychometrika, 76,} 306-317.

Cheung, M. W.-L. (2009). Constructing approximate confidence intervals for
parameters with structural constructing approximate confidence intervals for
parameters with structural equation models. \emph{Structural Equation Modeling,
16,} 267-294.

Pornprasertmanit, S., & Schneider, W. J. (2010). \emph{Efficient sample size for power and desired accuracy in Cohen's d estimation in two-group cluster randomized design} (Master Thesis). Illinois State University, Normal, IL. 

Pornprasertmanic, S., & Schneider, W. J. (2014). Accuracy in parameter estimation in cluster randomized designs. \emph{Psychological Methods}, \emph{19}, 356--379.}
\author{
Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
\dontrun{
# Examples for each function
ss.aipe.crd.es.nclus.fixedwidth(width=0.3, nindiv=20, es=0.5, estype=1, iccy=0.25, prtreat=0.5, 
	nrep=20)

ss.aipe.crd.es.nindiv.fixedwidth(width=0.3, 250, es=0.5, estype=1, iccy=0.25, prtreat=0.5, 
	nrep=20)

ss.aipe.crd.es.nclus.fixedbudget(budget=1000, nindiv=20, cluscost=0, indivcost=1, nrep=20, 
	prtreat=0.5, iccy=0.25, es=0.5)

ss.aipe.crd.es.nindiv.fixedbudget(budget=1000, nclus=200, cluscost=0, indivcost=1, nrep=20, 
	prtreat=0.5, iccy=0.25, es=0.5)

ss.aipe.crd.es.both.fixedbudget(budget=1000, cluscost=5, indivcost=1, es=0.5, estype=1, 
	iccy=0.25, prtreat=0.5, nrep=20)

ss.aipe.crd.es.both.fixedwidth(width=0.5, cluscost=5, indivcost=1, es=0.5, estype=1, iccy=0.25, 
	prtreat=0.5, nrep=20)

# Examples for different cluster size
ss.aipe.crd.es.nclus.fixedwidth(width=0.3, nindiv=20, es=0.5, estype=1, iccy=0.25, prtreat=0.5, 
nrep=20, diffsize = c(-2, 1, 0, 2, -1, 3, -3, 0, 0))

ss.aipe.crd.es.nclus.fixedwidth(width=0.3, nindiv=20, es=0.5, estype=1, iccy=0.25, prtreat=0.5, 
nrep=20, diffsize = c(0.6, 1.2, 0.8, 1.4, 1, 1, 1.1, 0.9))
}
}

