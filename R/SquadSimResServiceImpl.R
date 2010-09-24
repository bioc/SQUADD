# Class to create a table to be used as input for the prediction Map and the PCA
# 
# Author: admin
###############################################################################
library("RColorBrewer")

setClass(
		Class="SquadSimResServiceImpl",
		representation=representation(
				folder="character",
				time="numeric",
				ncolor="numeric",
				legend="vector",
				indexDeno="numeric",
				method="character",
				selectNode="character",
				conditionList="character"
		)
)


# constructor

simResService <- function(folder, time, ncolor=5, legend, indexDeno=1, method="lowess",selectNode="NA", conditionList="NA"){
	new(Class = "SquadSimResServiceImpl", ncolor=ncolor, folder=folder, time=time, legend=legend, indexDeno=indexDeno, method=method,selectNode=selectNode, conditionList=conditionList)
}

## method


## display the matrix of fitted value at a user defined time point
setMethod(
		f="getFittedTable",
		signature=signature(object="SquadSimResServiceImpl"),
		
		definition=function(object){
			listDir <- dir(object@folder,pattern = "*.txt", full.names=TRUE)
			vect <- lapply(listDir, getValues, time=object@time, meth=object@method)
			nm <- dir(object@folder, pattern = "*.txt",full.names=FALSE)
			mat <- as.data.frame(do.call("cbind", vect))
			names(mat) <- nm
			row.names(mat) <- object@legend
			return(mat)
		}
)

## Display a matrix of the SQUAD simulation for the user-selected nodes
setMethod( f= plotSimMatrix, signature = "SquadSimResServiceImpl", 
		definition = function (object){
			if(object@selectNode[1]=="NA"){
				cat("No Node selected. Please select node to display in the <selectNode> field")
				return(NULL)
			}
			filesDir <- dir(object@folder, pattern=".txt", full.names=TRUE)			
			## pdf(outputfile, width=15,height=10)
			par(mfrow=c(length(filesDir), length(object@selectNode)))
			a <- lapply(filesDir, plotByTab, selection=object@selectNode, legend=object@legend)
			## dev.off()
			## cat(outputfile,"\n")
		})


## Display a prediction map of the node activation changes between two conditions

setMethod(f="plotPredMap",
		signature = "SquadSimResServiceImpl",
		definition = function (object){
			
			matdf <- getFittedTable(object)
			
			# compute ratio
			matRatio <- matrix(nrow=nrow(matdf), ncol=(ncol(matdf)-1))
			colV <- 1:ncol(matdf)
			colV <- colV[-object@indexDeno]
			colN <- NULL
			denoV <- as.numeric(matdf[,object@indexDeno])+1
			ii <- 1
			for (i in colV){
				colNi <- paste(names(matdf)[i], "vs" , names(matdf)[object@indexDeno])
				colN <- c(colN,colNi)
				matRatio[,ii] <- (as.numeric(matdf[,i])+1)/denoV
				ii <- ii+1
			}
			
			matRatio <- as.data.frame(matRatio)
			names(matRatio) <- colN
			row.names(matRatio) <- row.names(matdf)
			getColdMap(df=matRatio,  ncolor=object@ncolor)
			
		})


## Display the correlation circle
setMethod(f="plotCC", signature = "SquadSimResServiceImpl",
		definition = function (object){
			if(object@conditionList[1]=="NA"){
				cat("Empty <conditionList> field\n")
				return(NULL)
			}
			mat <- getFittedTable(object)
			names(mat) <- object@conditionList
			res <- prcomp(mat,scale = TRUE)
			tmp <- strsplit(object@folder, "/")
			CCplot(res, tmp[[1]][length(tmp[[1]])])
		})



## function to get the fitted matrix
getValues <- function(filename, time, meth){
	line <- time*3
	tab <- read.table(filename)
	tab <- tab[,-1]
	col <- 1:length(tab)
	
	if(meth=="lowess") {
		dfL <- data.frame(do.call("cbind", lapply(col,getLowessValues, tab=tab)))
	} else if (meth=="linear") {
		dfL <- data.frame(do.call("cbind", lapply(col,getLeastSqValues, tab=tab)))
	}else {
		cat("fitting methods are 'lowess' &'linear'")
		return(NULL)
	}
	vectAtTime <-as.numeric(dfL[line,])
	return (vectAtTime)
}

getLowessValues <- function (x, tab, fit){	
	rowN <- nrow(tab)
	loFity <- lowess(1:rowN, tab[,x], f = 0.24)$y
	return(loFity)
}

getLeastSqValues <- function (x,tab){
	rowN <- nrow(tab)
	v <- 1:rowN
	w <- tab[,x]
	lsFity <- lm(w~v)$fitted.values
	return(lsFity)
}





## plot simulation matrix functions

plotByTab <- function (filename, selection, legend){
	cat(filename,"\n")
	tabcurr <- read.table(filename)
	t <- tabcurr[,1]
	tabcurr <- tabcurr[,-1] 
	names(tabcurr) <- legend
	
	if (length(selection)==1){
		getSimPlot(1,tab=tabcurr, time=t )
	}else {
		tabcurr <-  tabcurr[,selection]
		lapply(1:ncol(tabcurr), getSimPlot, tab= tabcurr, time=t)
	}
	
}



getSimPlot <- function(i, tab,time){
	
	y <- as.numeric(tab[,i])
	x <- time
	## checkEquals(length(x), length(y))
	colcurr <- palette()[i]
	
	lmFit <- lm(y~x)
	if(i == 1) {
		par(mar=c(2,3,1,1),font=2)
	}else {par(mar=c(2,1,1,1),font=2)}
	if(i==1){
		plot(x,y,lty=1, type = "b", pch = 1,ylim=c(0,1), axes=FALSE, col = colcurr)
		legend("right", legend=paste("y=", round(lmFit$coefficients[[2]], digits=3) ,"x +",round(lmFit$coefficients[[1]], digits=3)))
		
		axis(side=2, labels=NULL,tck=1, col ="light grey", cex.axis=1.5,font.axis=2)
		axis(side=1,at=NULL,labels=NULL, cex.axis=1.5,font.axis=2)
		abline(lmFit, col="red", lwd = 2)	
	}else {
		
		plot(x,y,lty=1, type = "b", pch = 1,ylim=c(0,1), axes=FALSE, col = colcurr)
		legend("right", legend=paste("y=", round(lmFit$coefficients[[2]], digits=3) ,"x +",round(lmFit$coefficients[[1]],  digits=3)))
		
		axis(side=2, at=NULL,labels=FALSE,tck=1, col ="light grey", cex.axis=1.5,font.axis=2)
		axis(side=1,at=NULL,labels=NULL, cex.axis=1.5,font.axis=2)
		abline(lmFit, col="red", lwd = 2)
		
	}
}









## plot prediction map

getColdMap <- function(df, outPicture=NULL, ncolor){
	df <- t(df)

	heatmap(data.matrix(df), Rowv=NA, Colv=NA, scale = "none", cexRow=1, cexCol=1,col=brewer.pal(ncolor, "Blues")
			,margins=c(10,15), na.rm=TRUE)
	
}



## plot correlation circle
CCplot <- function(x,title){
	xlabelpc1 <- paste("comp 1 (", round(summary(x)$importance[2,1]*100,1),"%)",sep="")
	ylabelpc2 <- paste("comp 2 (", round(summary(x)$importance[2,2]*100,1),"%)",sep="")
	a <- seq(0,2*pi,length=100)
	plot( cos(a), sin(a), type = 'l', lty=3, lwd=1,xlab = xlabelpc1, ylab = ylabelpc2, main = title)
	abline(h=0,v=0, lwd=1, lty=3)	
	n <- 2 #number PC
	df <- t(x$rotation)[1:n,]
	arrows(0,0,df[1,],df[2,],length=0.1)
	text(df[1,]-0.1,df[2,]-0.1,colnames(df),lwd=2)
}











## getter/setter

setMethod(f="[",
		signature="SquadSimResServiceImpl",
		definition=function(x,i,j,drop){
			switch(EXPR=i,
					"folder" = {return (x@folder)},
					"time"={return (x@time)},
					"ncolor"= {return (x@ncolor)}, "legend"={return(x@legend)}, 
					"indexDeno"={return(x@indexDeno)}, "method"={return(x@method)},
					"selectNode"={return(x@selectNode)},"conditionList"={return(x@conditionList)},
					stop("attribut does not exist"))
			
		}

)


setReplaceMethod(f="[",
		signature="SquadSimResServiceImpl",
		definition=function(x,i,j,value){
			switch(EXPR=i,
					"folder" = {x@folder <- value},
					"time"={x@time <- value}, 
					"ncolor"= {x@ncolor <- value}, "legend"={x@legend <- value}, 
					"indexDeno"={x@indexDeno <- value}, "method"={x@method <- value},
					"selectNode"={x@selectNode <- value},"conditionList"={x@conditionList <- value},
					stop("attribut does not exist"))
			
			validObject(x)
			return(x)
		}
)
























