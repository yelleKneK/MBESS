\name{HS}
\alias{HS}
\docType{data}
\title{ Complete Data Set of Holzinger and Swineford's (1939) Study}
\description{
The \emph{complete} data set of scores of 301 participants in 26 tests in Holzinger and Swineford's (1939) study.
}

\usage{data(HS)}
\format{
  A data frame with 301 observations on the following 34 variables.
  \describe{
    \item{\code{id}}{case number of participants (note there are skips) }
    \item{\code{sex}}{sex of participants}
    \item{\code{grade}}{grade in school of the participants with levels \code{Female} \code{Male}}
    \item{\code{age}}{the age (ignoring months into the year) of the participants}
    \item{\code{month_since_birthday}}{the number of months since the last birthday}
    \item{\code{age_months}}{age in months}
    \item{\code{age_years}}{age in years and months combined (more fine grained measure of years)}
    \item{\code{school}}{the school the participant is from with levels \code{Grant-White} \code{Pasteur}}
    \item{\code{t1_visual_perception}}{scores on visual perception test, test 1}
    \item{\code{t2_cubes}}{scores on cubes test, test 2}
    \item{\code{t3_paper_form_board}}{scores on paper form board test, test 3}
    \item{\code{t4_lozenges}}{scores on lozenges test, test 4}
    \item{\code{t5_general_information}}{scores on general information test, test 5}
    \item{\code{t6_paragraph_comprehension}}{scores on paragraph comprehension test, test 6}
    \item{\code{t7_sentence}}{scores on sentence completion test, test 7}
    \item{\code{t8_word_classification}}{scores on word classification test, test 8}
    \item{\code{t9_word_meaning}}{scores on word meaning test, test 9}
    \item{\code{t10_addition}}{scores on add test, test 10}
    \item{\code{t11_code}}{scores on code test, test 11}
    \item{\code{t12_counting_groups_of_dots}}{scores on counting groups of dots test, test 12}
    \item{\code{t13_straight_and_curved_capitals}}{scores on straight and curved capitals test, test 13}
    \item{\code{t14_word_recognition}}{scores on word recognition test, test 14}
    \item{\code{t15_number_recognition}}{scores on number recognition test, test 15}
    \item{\code{t16_figure_recognition}}{scores on figure recognition test, test 16}
    \item{\code{t17_object_number}}{scores on object-number test, test 17}
    \item{\code{t18_number_figure}}{scores on number-figure test, test 18}
    \item{\code{t19_figure_word}}{scores on figure-word test, test 19}
    \item{\code{t20_deduction}}{scores on deduction test, test 20}
    \item{\code{t21_numerical_puzzles}}{scores on numerical puzzles test, test 21}
    \item{\code{t22_problem_reasoning}}{scores on problem reasoning test, test 22}
    \item{\code{t23_series_completion}}{scores on series completion test, test 23}
    \item{\code{t24_woody_mccall}}{scores on Woody-McCall mixed fundamentals, form I test, test 24}
    \item{\code{t25_paper_form_board_r}}{scores on additional paper form board test, test 25}
    \item{\code{t26_flags}}{scores on flags test, test 26}
  }
}

\details{
Holzinger and Swineford (1939) data is widely cited, but generally only the Grant-White School data 
is used. The present dataset contains the complete data of Holzinger and Swineford (1939). 
 
A total number of 301 pupils, coming from Paster School and Grant-White School, who participated in Holzinger
and Swineford's (1939) study. This study consists of 26 tests, which are used to measure 
the participants' spatial, verbal, mental speed, memory, and mathematical ability. 

The spatial tests consist of \code{t1_visual_perception}, \code{t2_cubes}, \code{t3_paper_form_board}, \code{t4_lozenges}. Additional spatial tests are \code{t25_paper_form_board_r} (revised test 3) and \code{t26_flags}. \code{t25_paper_form_board_r} can (potentially) be used as a substitute for \code{t3_paper_form_board}. \code{t26_flags} is thought to be a possible substitute for \code{t4_lozenges}.
 
The verbal tests consist of \code{t5_general_information}, \code{t6_paragraph_comprehension}, \code{t7_sentence}, \code{t8_word_classification}, and \code{t9_word_meaning}.

The speed tests consist of \code{t10_addition}, \code{t11_code}, \code{t12_counting_groups_of_dots}, and \code{t13_straight_and_curved_capitals}.

The memory tests consist of \code{t14_word_recognition}, \code{t15_number_recognition}, \code{t16_figure_recognition}, \code{t17_object_number}, \code{t18_number_figure}, and \code{t19_figure_word}.

The mathematical-ability tests consist of \code{t20_deduction}, \code{t21_numerical_puzzles}, \code{t22_problem_reasoning}, 
\code{t23_series_completion}, and \code{t24_woody_mccall}.


}

\source{
Holzinger, K. J. and Swineford, F. A. (1939). A study in factor
analysis: The stability of a bi-factor solution. \emph{Supplementary Education
Monographs, 48}. University of Chicago.
}

\references{
Holzinger, K. J. and Swineford, F. A. (1939). A study in factor
analysis: The stability of a bi-factor solution. \emph{Supplementary Education
Monographs, 48}. University of Chicago.
}

\examples{
data(HS)
summary(HS)
}
\keyword{datasets}
