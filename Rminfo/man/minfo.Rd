\name{minfo}
\alias{minfo}
\title{Mutual Information Statistics}
\description{
    minfo calculates the percentage of mutual information for each factor.  It uses
the Shannon Information algorithm.
}
\usage{
    minfo(svars, ivars, iterate=FALSE,progress=TRUE)
}
\arguments{
    \item{svars}{State Variables. The state vector.}
    \item{ivars}{Infomation Varibles. This is the data.  It must be descrete.
        If each column is not already a factor it will be converted using as.factor().  Take into
        account how you want NA's to be treated when creating the factors.}
    \item{iterate}{If true, minfo will calculate the mutual information for each variable
        then combine the most informative varible with all the others and recalculate
        the mutual information for all of the tuples.  This will continue with the "best"
        triple, quad, etc until some combination completely predicts the state vector.  This
        provides a "greedy" method of finding good combinations.}
    \item{progress}{If true progress is reported at each step.}
}
\value{
    If iterate=FALSE, the mutual information for each variable is return.  If TRUE, 
    the "best" vector of predictors is returned using the greedy algorithm described above.
    The names of the vector enumerate the columns of "ivars" while the values are the
    corresponding combination's mutual information.}

\author{ 
    James R. Oliphant \email{olijam@netzero.com}}
    
\seealso{gom}

\examples{
    ## The anorexia data set consists of a survey of 500 teen-age women.  
    ## Variable 95 (DIET14) is the question:
    ##  "I make myself throw up to help with weight control."
    ## Answers consist of: 0--never, 1--rarely , 2--sometimes, 3--often, 4--almost always.
    ##  See anorexia.README for the questions in the survey.
    data(anorexia)
    minfo(anorexia[,95],anorexia)

    ## "spore_color", "stalk_rings", "stalk_color" completely predicts the "genius" of mushrooms
    ## The tuple "spore_color","stalk_rings" specifies 93.6% of "genius"
    data(mushrooms)
    minfo(mushrooms$genius,mushrooms[,6:35],iterate=TRUE,progress=TRUE)
}

\keyword{minfo Information Shannon}
