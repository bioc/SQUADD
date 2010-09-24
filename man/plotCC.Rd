\name{plotCC}
\alias{plotCC}
\title{Generate the Correlation Circle.}
\description{
This method generates the correlation circle. Correlation circle is used to assess the coherence of a model. Indeed, Each vector represents our conditions and the angle between each vector give the correlation between the variables. if the angle is higher than 90 degrees there is no correlation.
}
\usage{
plotCC(object)
}
\arguments{
  \item{object}{
	Object of the Class \code{SquadSimResServiceImpl}.
}
}
\details{
The percentage of variability explained by the two first principal components are displayed in parenthesis in the axis label. 
}

\author{
Martial Sankar
}

\seealso{
More specific information in \code{\link{plotCC-methods}}.
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
		method="lowess",
		conditionList=c("Normal", "brxlof", "BRrescue","brxarfilof", 
		             "brxarfiBRrescue", "brxgof")
)
## call method
plotCC(obj)
}

\keyword{correlation}
\keyword{PCA}
\keyword{circle}
