\name{getFittedTable}
\alias{getFittedTable}
\title{Build the table of fitted values.}
\description{
This method creates a table containing the interpolated/fitted value of the simulation curve at the user-defined time point. Each column represents the condition and each line, the component node. It is used as basis for the Correlation Circle and the Prediction Map.
}
\usage{
getFittedTable(object)
}
\arguments{
  \item{object}{
	Object of the Class \code{SquadSimResServiceImpl}.
}
}
\details{
This method firstly interpolate the SQUAD Simulation curve with a global (linear) or local (lowess) interpolation. A local interpolation is more sensitive to change in the curve when applying perturbation or pulse.
}
\value{
returns a matrix of size length of component x length of conditions. 
}
\references{
(1) Cleveland, W. S. (1979) Robust locally weighted regression and smoothing scatterplots. J. Amer. Statist. Assoc. 74, 829:836.
(2) Chambers, J. M. (1992) Linear models. Chapter 4 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth And Brooks/Cole. 
}
\author{
Martial Sankar
}

\seealso{
More specific information in \code{\link{getFittedTable-methods}}.
For more informations about the input object see \code{\link{SquadSimResServiceImpl-class}}.
}
\examples{
fpath <- system.file("extdata", package="SQUADD")
fileModel <- file.path(fpath,"data_IAA")
nm <- c("ARF(a)", "ARF(i)", "AR_Genes", "Aux/IAA", "BES1/BZR1",
 "BIN2", "BR", "BRI1-BAK1","BRR_Genes","BRX","BR_Biosynthesis","BZR1",
  "DO", "IAA", "IAA_Biosynthesis", "NGA1", "PIN", "SCFTir1",
  "StimAux", "StimBR")
t <- 50

## call constructor
obj <- simResService(
		folder=fileModel,
		time=t,
		ncolor=5,
		legend=nm,
		indexDeno=1,
		method="lowess")

## Apply method 
tab <- getFittedTable(obj)
}

\keyword{fit}
\keyword{interpolate}