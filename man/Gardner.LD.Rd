\name{Gardner.LD}
\alias{Gardner.LD}
\docType{data}
\title{The Gardner learning data, which was used by L.R. Tucker}
\description{
Repeated measures data on 24 participants, each with 21 trials (each trial based on 20 replications).
}
\usage{data(Gardner.LD)}
\format{
  A data frame where the rows represent the timepoints for the individuals.
  \describe{
    \item{\code{ID}}{: a numeric vector}
    \item{\code{Trial}}{: a numeric vector}
    \item{\code{Score}}{: a numeric vector}
    \item{\code{Group}}{: a numeric vector}
  }
}
\details{
The 24 participants of this study were presented with 420 presentations of four letters 
where the task was to identify the next letter that was to be presented. Twelve of the participants (Group 1) were
presented the letters S, L, N, and D with probabilities .70, .10, .10, and .10, respectively.
The other 12 participants (Group 2) were presented the letter L with probability .70 and three other
letters, each with a probability of .10. The 420 presentations were (arbitrarily it seems) 
grouped into 21 trials of 20 presentations. The score for each trial was the number of times the individual
correctly guessed the dominant letter. The participants were naive to the probability that the
letters would be presented. Other groups of individuals (although the data is not available) 
were tested under a different probability structure. The data given here is thus known as the 70-10-10-10 group from 
Gardner's paper. L. R. Tucker used this data set to illustrate methods for understanding change.
}

\source{Tucker, L. R. (1960). Determination of Generalized Learning Curves by Factor Analysis,
Educational Testing Services, Princeton, NJ.
}
\references{
Gardner, R. A., (1958). Multiple-choice decision-behavior, \emph{American Journal of Psychology}, 71,
710--717.
}
%\examples
%{
%}
\keyword{datasets}
