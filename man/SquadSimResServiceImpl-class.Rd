\name{SquadSimResServiceImpl-class}
\Rdversion{1.1}
\docType{class}
\alias{SquadSimResServiceImpl-class}
\alias{SquadSimResServiceImpl}
\alias{simResService}
\alias{[<-,SquadSimResServiceImpl-method}
\alias{[,SquadSimResServiceImpl-method}

\title{Class "SquadSimResServiceImpl"}
\description{
The Class \code{SquadSimResServiceImpl} contains the information needed to plot the simulation matrix, the prediction map and the correlation circle.
}

\section{Slots}{
  \describe{
    \item{\code{folder}:}{Object of class \code{"character"}, path to the folder to find the SQUAD (1) data.}
    \item{\code{time}:}{Object of class \code{"numeric"}, time point used to get the interpolated activate state value.}
    \item{\code{ncolor}:}{Object of class \code{"numeric"}, number of color code to used in the prediction map. default = 5 }
    \item{\code{legend}:}{Object of class \code{"vector"}, contains the name of each node of the model (Must appear in the same order than the simulation). }
    \item{\code{indexDeno}:}{Object of class \code{"numeric"}, the index of the file to used as denominator (the names of each file in the folder should contain an index, i.e 1_xxx.txt), default = 1 }
    \item{\code{method}:}{Object of class \code{"character"}, the fitting method to interpolate the curve ("linear" to get a straight line (2) or "lowess" to get a smoothed interpolation (3)).}
    \item{\code{selectNode}:}{Object of class \code{"character"}, contains the name of the node to display.}
    \item{\code{conditionList}:}{Object of class \code{"character"}, list of the perturbation condition names.}
    }
    
}

\section{Constructor}{
\code{simResService(folder, time, ncolor, legend, indexDeno, method, selectNode, conditionList)}: construct a instance of \code{SquadSimResServiceImpl}.
}

\section{Methods}{
  \describe{
    \item{getFittedTable}{\code{signature(object = "SquadSimResServiceImpl")}: Build and display the table of interpolated/fitted values (see \link{getFittedTable} for more details).}
    \item{plotSimMatrix}{\code{signature(object = "SquadSimResServiceImpl")}: plot the simulation matrix (see \link{plotSimMatrix} for more details).}
    \item{plotPredMap}{\code{signature(object = "SquadSimResServiceImpl")}: plot the prediction map (see \link{plotPredMap} for more details).}
    \item{plotCC}{\code{signature(object = "SquadSimResServiceImpl")}: plot the correlation circle from PCA analyis (see \link{plotCC} for more details).}
	 }
}

\section{Accessor methods}{
See examples Section.
}

\references{
(1) Di Cara A, Garg A, De Micheli G, Xenarios I, Mendoza L (2007) Dynamic simulation of regulatory networks using SQUAD. BMC Bioinformatics 8: 462.
(2) Cleveland, W. S. (1979) Robust locally weighted regression and smoothing scatterplots. J. Amer. Statist. Assoc. 74, 829:836.
(3) Chambers, J. M. (1992) Linear models. Chapter 4 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth And Brooks/Cole. 
}
\author{
Martial Sankar
}

\examples{
fpath <- system.file("extdata", package="SQUADD")
fileModel <- file.path(fpath,"data_IAA")
nm <- c("ARF(a)", "ARF(i)", "AR_Genes", "Aux/IAA", "BES1/BZR1", "BIN2", "BR", "BRI1-BAK1","BRR_Genes","BRX","BR_Biosynthesis","BZR1", "DO", "IAA", "IAA_Biosynthesis", "NGA1", "PIN", "SCFTir1","StimAux", "StimBR")
t <- 50

## call constructor
obj <- simResService(folder=fileModel, time=t, ncolor=5,legend=nm, indexDeno=1,	method="lowess")

## Apply method 
tab <- getFittedTable(obj)
obj["selectNode"] <-c("DO","IAA_Biosynthesis","BR_Biosynthesis", "IAA", "BR")

plotSimMatrix(obj)
plotPredMap(obj)

obj["conditionList"] <- c("Normal", "brxlof", "BRrescue","brxarfilof", "brxarfiBRrescue", "brxgof")
plotCC(obj)

 ## getter/setter test
obj["time"] <- 45
}
\keyword{classes}
\keyword{constructor}
\keyword{accessor}
\keyword{setter}
\keyword{getter}