\name{intr.plot.2d}
\alias{intr.plot.2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting Conditional Regression Lines with Interactions in Two Dimensions }
\description{
To plot regression lines for one two-way interactions, holding one of the predictors (in this function, z) at 
values -2, -1, 0, 1, and 2 standard deviations above the mean.
}
\usage{

intr.plot.2d(b.0, b.x, b.z, b.xz,x.min=NULL, x.max=NULL, x=NULL, 
n.x=50, mean.z=NULL, sd.z=NULL, z=NULL,xlab="Value of X",  
ylab="Dependent Variable", sd.plot=TRUE, sd2.plot=TRUE, sd_1.plot=TRUE, 
sd_2.plot=TRUE, type.sd=2, type.sd2=3, type.sd_1=4, type.sd_2=5, 
legend.pos="bottomright", legend.on=TRUE, ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b.0}{ the intercept }
  \item{b.x}{ regression coefficient for predictor x }
  \item{b.z}{ regression coefficient for predictor z  }
  \item{b.xz}{ regression coefficient for the interaction of predictors x and z }
  \item{x.min,x.max}{the range of x used in the plot}
  \item{x}{ a specific predictor vector x, used instead of \code{x.min} and \code{x.max} }
  \item{n.x}{number of elements in predictor vector x}
  \item{mean.z}{ mean of predictor z }
  \item{sd.z}{ standard deviation of predictor z }
  \item{z}{a specific predictor vector z, used instead of \code{z.min} and \code{z.max} }
  \item{xlab}{ title for the axis which the predictor x is on }
  \item{ylab}{ title for the axis which the dependent y is on }
  \item{sd.plot, sd2.plot, sd_1.plot, sd_2.plot}{ whether or not to plot 
  the regression line holding z at values 1, 2, -1, and -2 standard deviations above the mean, respectively. 
  Default values are all \code{TRUE}. }
  \item{type.sd, type.sd2, type.sd_1, type.sd_2}{ types of lines to be plotted holding z at values 1, 2,
  -1, and -2 standard deviations above the mean, respectively. Default are line type 2,3,4, and 5, respectively. }
  \item{legend.pos}{ position of the legend; possible options are \code{"bottomright"}, \code{"bottom"}, 
  \code{"bottomleft"}, \code{"left"}, \code{"center"}, \code{"right"}, \code{"topleft"}, 
  \code{"top"}, and \code{"topright"}. }
  \item{legend.on}{whether or not to show the legend}
  \item{\dots}{ allows one to potentially include parameter values for inner functions }

}
\details{
To input the predictor x, one can use either the limits of x (\code{x.max} and \code{x.min}) , or a specific vector x (\code{x}). 
To input the predictor z, one can use either the mean and standard deviation of z (\code{mean.z} and \code{sd.z} ), or a specific vector z (\code{z}).   
}
%\value{}
\references{Cohen, J., Cohen, P., West, S. G. and Aiken, L. S. (2003). \emph{Applied multiple regression/correlation 
analysis for the behavioral sciences} (3rd ed.). Mahwah, NJ: Erlbaum.}


\author{ Keke Lai, Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})  }

\note{ Sometimes some of the regression lines are outside the default scope of the coordinates 
and thus cannot be seen; in such situations, one needs to, by entering additional arguments, adjust the scope to
let proper sections of regression lines be seen. Refer to examples below for more details.
}
\seealso{ \code{intr.plot} }
\examples{
## A situation where one regression line is outside the default scope of the coordinates
intr.plot.2d(b.0=16, b.x=2.2, b.z=2.6, b.xz=.4, x.min=0, x.max=20, mean.z=0, sd.z=3)

## Adjust the scope of x and y axes so that proper sections of regression lines can be seen 
intr.plot.2d(b.0=16, b.x=2.2, b.z=2.6, b.xz=.4, x.min=0, x.max=50, mean.z=0, 
sd.z=3, xlim=c(0,50), ylim=c(-20,100) )

## Use specific vector(s) to define the predictor(s) 
intr.plot.2d(b.0=16, b.x=2.2, b.z=2.6, b.xz=.4, x=c(1:10), z=c(0,2,4,6,8,10))

intr.plot.2d(b.0=16, b.x=2.2, b.z=2.6, b.xz=.4, x.min=0, x.max=20, 
z=c(1,3,6,7,9,13,16,20), ylim=c(0,100))

## Change the position of the legend so that it does not block regression lines
intr.plot.2d(b.0=10, b.x=-.3, b.z=1, b.xz=.5, x.min=0, x.max=40, mean.z=-5, sd.z=3, 
ylim=c(-100,100),legend.pos="topright" )

}
\keyword{ regression }% at least one, from doc/KEYWORDS
