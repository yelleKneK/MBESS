\name{mediation.effect.plot}
\alias{mediation.effect.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Visualizing mediation effects}
\description{ Create a mediation effect plot
}
\usage{
mediation.effect.plot(x, mediator, dv, ylab = "Dependent Variable", 
xlab = "Mediator", main = "Mediation Effect Plot", 
pct.from.top.a = 0.05, pct.from.left.c = 0.05, arrow.length.a = 0.05, 
arrow.length.c = 0.05, legend.loc = "topleft", file = "", pch = 20, 
xlim = NULL, ylim = NULL, save.pdf = FALSE, save.eps = FALSE, 
save.jpg = FALSE, ...)
}

\arguments{
\item{x}{vector of the predictor/independent variable}
\item{mediator}{vector of the mediator variable}
\item{dv}{vector of the dependent/outcome variable}
  \item{ylab}{ \emph{y}-axis title label }
  \item{xlab}{ \emph{x}-axis title label }
  \item{main}{ main title label }
  \item{pct.from.top.a}{figure fine tuning adjustment}
  \item{pct.from.left.c}{figure fine tuning adjustment}
  \item{arrow.length.a}{figure fine tuning adjustment }
  \item{arrow.length.c}{ figure fine tuning adjustment}
  \item{legend.loc}{ specify the location of the legend }
  \item{file}{ file name of the plot to be saved (not necessary) }
  \item{pch}{ plotting character }
  \item{xlim}{ limits for the \emph{x}-axis }
  \item{ylim}{ limits for the \emph{y}-axis }
  \item{save.pdf}{\code{TRUE} or \code{FALSE} if the produced figure should be saved as a PDF file}
  \item{save.eps}{\code{TRUE} or \code{FALSE} if the produced figure should be saved as an EPS file}
  \item{save.jpg}{\code{TRUE} or \code{FALSE} if the produced figure should be saved as a JPG file}
  \item{\dots}{to incorporate options from interval functions}
}
\details{


Merrill (1994; see also MacKinnon, 2008; MacKinnon et al., 2007; Sy, 2004) presents a method that involves plotting the indirect effect as the vertical distance between two lines. Fritz and MacKinnon (2008) present a detailed exposition of this method too. Preacher and Kelley (2011) discuss this plotting method and implement their own code, which was also independently done as part of Fritz and MacKinnon (2008).

In this type of plot, the two horizontal lines correspond to the predicted values of Y regressed on X at the mean of X and at one unit above the mean of X. The distance between these two lines is thus \eqn{\hat{c}}. The two vertical lines correspond to predicted values of M regressed on X at the same two values of X. The distance between these lines is \eqn{\hat{a}}. The lines corresponding to the regression of Y on M (controlling for X) are plotted for the same two values of X. 

}
\value{
A figure is returned.
}
\references{
Fritz, M. S., & MacKinnon, D. P. (2008). A graphical representation of the mediated effect. \emph{Behavior Research Methods, 40}, 55--60.

MacKinnon, D. P. (2008). \emph{Introduction to statistical mediation analysis}. Mahwah, NJ: Erlbaum.

MacKinnon, D. P., Fairchild, A. J., & Fritz, M. S. (2007). Mediation analysis. \emph{Annual Review of Psychology, 58}, 593--614.

Merrill, R. M. (1994). \emph{Treatment effect evaluation in non-additive mediation models}. Unpublished dissertation, Arizona State University.

Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation models: Quantitative and graphical strategies for communicating indirect effects. \emph{Psychological Methods}, \emph{16}, 93--115.

Sy, O. S. (2004). \emph{Multilevel mediation analysis: Estimation and applications}. Unpublished dissertation, Kansas State University.
}
\author{Ken Kelley (University of Notre Dame; KKelley@nd.edu)}

\note{Requires raw data.}

\seealso{\code{\link{mediation.effect.plot}}, \code{\link{mediation.effect.bar.plot}}}
%\examples{}

\keyword{multivariate}