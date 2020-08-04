\name{intr.plot}
\alias{intr.plot}

\title{  Regression Surface Containing Interaction }
\description{
To plot a three dimentional figure of a multiple regression surface containing one two-way interaction.
}
\usage{
intr.plot(b.0, b.x, b.z, b.xz, x.min = NULL, x.max = NULL, z.min = NULL, 
z.max = NULL, n.x = 50, n.z = 50, x = NULL, z = NULL, col = "lightblue", 
hor.angle = -60, vert.angle = 15, xlab = "Value of X", zlab = "Value of Z", 
ylab = "Dependent Variable", expand = 0.5, lines.plot=TRUE, col.line = "red", 
line.wd = 2, gray.scale = FALSE, ticktype="detailed", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{b.0}{ the intercept }
  \item{b.x}{ regression coefficient for predictor x }
  \item{b.z}{ regression coefficient for predictor z }
  \item{b.xz}{ regression coefficient for the interaction of predictors x and z }
 \item{x.min, x.max, z.min, z.max}{ranges of x and z. The regression surface defined by these limits will be plotted.}
\item{n.x}{ number of elements in predictor vector x; number of points to be plotted on the regression surface; default is 50}
  \item{n.z}{ number of elements in predictor vector z; number of points to be plotted on the regression surface; default is 50}
 \item{x}{ a specific predictor vector \code{x}, used instead of \code{x.max} and \code{x.min}}
  \item{z}{ a specific predictor vector \code{z}, used instead of \code{z.max} and \code{z.min} }
   \item{col}{ color of the regression surface; default is lightbule}
 \item{hor.angle}{ rotate the regression surface horizontally; default is -60 degree }
  \item{vert.angle}{ rotate the regression surface vertically; default is 15 degree }
 \item{xlab}{ title for the axis which the predictor \code{x} is on  }
  \item{zlab}{ title for the axis which the predictor \code{z} is on }
  \item{ylab}{ title for the axis which the dependent \code{y} is on }
  \item{expand}{ default is 0.5; expansion factor applied to the axis of the dependent variable. Often used with 0 < \code{expand} < 1 to shrink the plotting box in the direction of the dependent variable's axis. }
  \item{lines.plot}{ whether or not to plot on the regression surface regression lines holding z at
  values 0, 1, -1, 2, -2 above the mean; default is \code{TRUE}. }
  \item{col.line}{ the color of regression lines plotted on the regression surface; default is red }
  \item{line.wd}{ the width of regression lines plotted on the regression surface; default is 2 }
  \item{gray.scale}{ whether or not to plot the figure black and white; default is \code{FALSE} }
 \item{ticktype}{ whether the axes should be plotted with (\code{"detailed"}) or without (\code{"simple"}) tick marks}
 \item{\dots}{ allows one to potentially include parameter values for inner functions }
}

\details{
The user can input either the limits of \code{x} and \code{z}, or specific \code{x} and \code{z} vectors, to draw the regression surface. If the user inputs simply the limits of the predictors, the function would generate predictor vectors for plotting. 
If the user inputs specific predictor vectors, the function would plot the regression surface based on those vectors. 
}

%\value{}

\references{Cohen, J., Cohen, P., West, S. G. and Aiken, L. S. (2003). \emph{Applied multiple regression/correlation analysis for the behavioral sciences} (3rd ed.). Mahwah, NJ: Erlbaum. }

\author{Keke Lai (University of California -- Merced) and Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }

\note{ If the user enters specific vectors instead of the ranges of predictors, please make sure
elements in those vectors are in ascending order. This is required by function \code{persp}, which 
is used within this function.  
}
\seealso{ \code{intr.plot.2d}, \code{persp} }
\examples{
##  A way to replicate the example given by Cohen et al. (2003) (pp. 258--263):
##  The regression equation with interaction is y=.2X+.6Z+.4XZ+2
##  To plot a regression surface and regression lines of Y on X holding Z 
##  at -1, 0, and 1 standard deviation above the mean

x<- c(0,2,4,6,8,10)
z<-c(0,2,4,6,8,10)
intr.plot(b.0=2, b.x=.2, b.z=.6, b.xz=.4, x=x, z=z)

## input limits of the predictors instead of specific x and z predictor vectors
intr.plot(b.0=2, b.x=.2, b.z=.6, b.xz=.4, x.min=5, x.max=10, z.min=0, z.max=20)

intr.plot(b.0=2, b.x=.2, b.z=.6, b.xz=.4, x.min=0, x.max=10, z.min=0, z.max=10, 
col="gray", hor.angle=-65, vert.angle=10)

##  To plot a black-and-white figure
intr.plot(b.0=2, b.x=.2, b.z=.6, b.xz=.4, x.min=0, x.max=10, z.min=0, z.max=10, 
gray.scale=TRUE)

## to adjust the tick marks on the axes
intr.plot(b.0=2, b.x=.2, b.z=.6, b.xz=.4, x.min=0, x.max=10, z.min=0, z.max=10, 
ticktype="detailed", nticks=8)
}


\keyword{ regression }
