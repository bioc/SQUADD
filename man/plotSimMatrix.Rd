\name{plotSimMatrix}
\alias{plotSimMatrix}
\title{Generate the SQUAD Simulation matrix display.}
\description{
This method is used to create a simulation matrix. Each columns correspond to the user-selected component node and the each rows, to the conditions.
}
\usage{
plotSimMatrix(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of the Class \code{SquadSimResServiceImpl}.
}
}

\author{
Martial Sankar
}

\seealso{
More specific information in \code{\link{plotSimMatrix-methods}}.
For more informations about the input object see \code{\link{SquadSimResServiceImpl-class}}.
}
\examples{
fpath <- system.file("extdata", package="SQUADD")
fileModel <- file.path(fpath,"data_IAA")
nm <- c("ARF(a)", "ARF(i)", "AR_Genes", "Aux/IAA", "BES1/BZR1",
 "BIN2", "BR", "BRI1-BAK1","BRR_Genes","BRX","BR_Biosynthesis","BZR1",
  "DO", "IAA", "IAA_Biosynthesis", "NGA1", "PIN", "SCFTir1","StimAux", "StimBR")
t <- 50

## call constructor
obj <- simResService(
		folder=fileModel,
		time=t,
		ncolor=5,
		legend=nm,
		indexDeno=1,
		method="lowess",
		selectNode=c("DO","IAA_Biosynthesis","BR_Biosynthesis", "IAA", "BR"))

## call method
plotSimMatrix(obj)		
		
		
		
}
\keyword{simulation}

