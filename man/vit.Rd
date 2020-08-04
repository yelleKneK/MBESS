\name{vit}
\alias{vit}
\title{ Visualize individual trajectories  }
\description{
  A function to help visualize individual trajectories in a longitudinal (i.e., analysis of change) context.
}
\usage{
vit(id = "", occasion = "", score = "", Data = NULL, group = NULL, 
subset.ids = NULL, pct.rand = NULL, number.rand = NULL, 
All.in.One = TRUE, ylab = NULL, xlab = NULL, same.scales = TRUE, 
plot.points = TRUE, save.pdf = FALSE, save.eps = FALSE,
 save.jpg = FALSE, file = "", layout = c(3, 3), col = NULL, 
 pch = 16, cex = 0.7, ...)
}
\arguments{
  \item{id}{ string variable of the column name of id }
  \item{occasion}{ string variable of the column name of time variable }
  \item{score}{ string variable of the column name where the score (i.e., dependent variable) is located }
  \item{Data}{ data set with named column variables (see above) }
  \item{group}{ if plotting parameters should be conditional on group membership }
  \item{subset.ids}{ id values for a selected subset of individuals }
  \item{pct.rand}{ percentage of random trajectories to be plotted }
  \item{number.rand}{ number of random trajectories to be plotted }
  \item{All.in.One}{ should trajectories be in a single or multiple plots }
  \item{ylab}{ label for the ordinate (i.e., y-axis; see par) }
  \item{xlab}{ label for the abscissa (i.e., x-axis; see par) }
  \item{same.scales}{ should the y-axes have the same scales }
  \item{plot.points}{ should the points be plotted }
  \item{save.pdf}{ save a pdf file }
  \item{save.eps}{ save a postscript file }
  \item{save.jpg}{ save a jpg file }
  \item{file}{ file name and file path for the graph(s) to save, if \code{file=""} a file would be saved in the current working directory }
  \item{layout}{ define the per-page layout when \code{All.in.One=FALSE} }
  \item{col}{ color(s) of the line(s) and points }
  \item{pch}{ plotting character(s); see par }
  \item{cex}{ size of the points (1 is the R default; see par) }
  \item{\dots}{ optional plotting specifications }
}
\details{
  This function makes visualizing individual trajectories simple. 
  Data should be in the "univariate format" (i.e., the same format as lmer and nlme data.)
}
\value{
  Returns a plot of individual trajectories with the specifications provided.
}

%\references{ }

\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) and Po-Ju Wu (Indiana University)}
%\note{ }

\seealso{ par, nlme, vit.fitted,  }

\examples{
\dontrun{
data(Gardner.LD)

# Although many options are possible, a simple call to
# 'vit' is of the form:
# vit(id="ID", occasion= "Trial", score= "Score", Data=Gardner.LD)

# Now color is conditional on group membership.
# vit(id="ID", occasion= "Trial", score="Score", Data=Gardner.LD, 
# group="Group")

# Now randomly selects 50% of the sample to plot
# vit(id="ID", occasion= "Trial", score="Score", Data=Gardner.LD, 
# pct.rand=50, group="Group")

# Specified individuals are plotted (by group)
# vit(id="ID", occasion= "Trial", score="Score", Data=Gardner.LD, 
# subset.ids=c(1, 4, 8, 13, 17, 21), group="Group")

# Now colors for groups are changed .
# vit(id="ID", occasion= "Trial", score="Score", Data=Gardner.LD, 
# group="Group",subset.ids=c(1, 4, 8, 13, 17, 21), col=c("Green", "Blue"))

# Now each individual specified is plotted separately.
# vit(id="ID", occasion= "Trial", score="Score", Data=Gardner.LD, 
# group="Group",subset.ids=c(1, 4, 8, 13, 17, 21), col=c("Green", "Blue"),
# All.in.One=FALSE)
}
}

\keyword{hplot}
\keyword{device}
\keyword{dynamic}
