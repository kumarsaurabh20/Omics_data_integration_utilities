panelFunction <- function(x,y) { 
	require(lattice)
	panel.xyplot(x,y,cex=.2,pch=19, xlab="A",ylab="M")
	left <- -21.5
	right <- -16.25
	yamt <- 3
	m <- (right-left)/yamt
	b1 <- -m*left
	b2 <- -m*right
	mdpt <- (left+right)/2
	otherleft <- -23.9
	otherright <- -24.5
	panel.text(-17,1, round(100*mean( y < (m*x+b1) & y > (m*x+b2) & y < 0 & y > -yamt ),2),col="black",cex=1.5)
	panel.text(-22.5,4, round(mean(x<-22.9 & y>2)*100,2), col="red", cex=1.5)
	panel.text(-22.5,-4, round(100*mean(x<-22.9 & y< -2),2), col="orange",cex=1.5)
	panel.grid()
	panel.lines( c(right,mdpt,left,mdpt,right), c(0,-yamt,-yamt,0,0), lwd=3,col="black" ) 
	panel.lines( c(otherright,otherleft,otherleft,otherright), c(2,2,5.6,5.6,2), lwd=3,col="red" ) 
	panel.lines( c(otherright,otherleft,otherleft,otherright), c(-2,-2,-5.6,-5.6,-2), lwd=3,col="orange" ) 
}


makeROCCurve <- function(depth=30, low.diff=.1, high.diff=.3, xlim=c(0,.2), 
lwds=3, cols=c("blue","red","salmon","orange","darkgreen"), labels=c("Naive","ABCD-DNA (SNP 6.0)","ABCD-DNA (gDNA-seq)","RSEG","DiffBind"),
num.regions=50, min.cpgw=10, lcex=1, symmetric=FALSE) {
	
	stopifnot( length(cols)==length(labels) )
	stopifnot( length(labels)==length(scoresList) )
	
	nCurves <- length(cols)
	if( length(lwds) == 1 ) lwds <- rep(lwds, nCurves)
    stopifnot( length(lwds)==length(scoresList) )
	
	require(ROCR)
	tf <- rep(NA, nrow(betavals))
	md <- betavals[,2]-betavals[,1]
	tf[abs(md) > high.diff] <- 1
	tf[abs(md) < low.diff] <- 0
	if( symmetric ) {
		for(cn in cngrp) {
			pos <- !is.na(tf) & tf==1 & md > 0 & grp==cn
			neg <- !is.na(tf) & tf==1 & md < 0 & grp==cn
			cat(cn,"pos",sum(pos),"\n")
			cat(cn,"pos",sum(neg),"\n")
			pmn <-  sum(pos)-sum(neg)
			if( pmn >= 0 ) {
# remove some of the positive direction ones
                s <- sample( which(pos), pmn )
			} else {
# remove some of the negative direction ones
                s <- sample( which(neg), -pmn )
			}
			tf[s] <- -99
		}
	}
	tf[!(grp %in% cngrp)] <- -99
	tf[cpgdens < min.cpgw] <- -99
	
print(table(tf,sign(md),grp))
	
# calculate ROC curves
	keep <- which(tf >= 0)
	
	scores <- lapply(scoresList, function(u) abs(u[keep]))
#scores <- list(naive=abs(glmz[keep]),adjusted=abs(glmadjz[keep]),gdna=abs(glmadjz1[keep]),rseg=sx[keep],db=abs(dbz[keep]))
	scores <- lapply(scores, function(u) split(u, grp[keep]))
	scores <- unlist(scores, recursive=FALSE)
	groups <- lapply(scoresList, function(u) tf[keep])
#groups <- list(naive=tf[keep],adjusted=tf[keep],gdna=tf[keep],sx=tf[keep],db=tf[keep])
	groups <- lapply(groups, function(u) split(u, grp[keep]))
	groups <- unlist(groups, recursive=FALSE)
	
	pred <- prediction(scores, groups)
	perf <- performance(pred,"tpr", "fpr")
	
    nm <- names(scores)
	
	names(perf@alpha.values) <- names(perf@x.values) <- names(perf@y.values) <- nm
	n <- length(scores)
	
	lwds <- rep(lwds,each=n/nCurves)
	cols <- rep(cols,each=n/nCurves)	
	print(cols)
	slen <- sapply(scores,length)
	
	for(i in 1:(n/nCurves)) {
		ind <- c(i+(0:(nCurves-1))*n/nCurves)
		mn <- strsplit(names(scores)[ind],"\\.")
		idx <- .getSmallSeq(perf@x.values[[i]])
		plot( perf@x.values[[i]][idx], perf@y.values[[i]][idx], lwd=lwds[i], type="l", cex.main=1.4, cex.lab=1.4, cex.axis=1.2,
			 col=cols[i], xlab="FPR", ylab="TPR", xlim=xlim, main=paste(mn[[1]][2]," (",slen[i],")",sep=""))
		for(j in 0:(nCurves-1)) {
			this <- i+j*n/nCurves
			idx <- .getSmallSeq(perf@x.values[[this]])
			lines( perf@x.values[[this]][idx], perf@y.values[[this]][idx],  lwd=lwds[this], col=cols[this])
		}		
		legend("bottomright", labels, lwd=lwds[ind], col=cols[ind], cex=lcex)
	}
	a <- performance(pred,"auc")@y.values
	str(a)
	invisible(list(tf=tf,pred=pred,perf=perf,scores=scores,groups=groups))
}


.getSmallSeq <- function(xvals, max.length=200) {
  npoints <- length(xvals)
  if (npoints > max.length)
    idx <- seq(1,npoints,length=max.length)
  else
    idx <- 1:npoints
}




makeROCCurve.overall <- function(main,md,xlim=c(0,1),symmetric=TRUE,high.diff=0.4, low.diff=0.1, 
cols=c("blue","darkgreen","red"), labels=c("Naive","BATMAN","ABCD-DNA")) {
	
	require(ROCR)
	tf <- rep(-99, length(m27kGR))
	
	tf[abs(md) > high.diff] <- 1
	tf[abs(md) < low.diff] <- 0
	if( symmetric ) {
		pos <- !is.na(tf) & tf==1 & md > 0
		neg <- !is.na(tf) & tf==1 & md < 0
		pmn <-  sum(pos)-sum(neg)
		
		if( pmn >= 0 ) {
			s <- sample( which(pos), pmn )
		} else {
			s <- sample( which(neg), -pmn )
		}
		tf[s] <- -99
	}

	print(table( tf, sign(md) ))

	keep <- which(tf >= 0)
	
	scores <- lapply(scoresList, function(u) abs(u[keep]))
	groups <- lapply(scoresList, function(u) tf[keep])
	
	pred <- prediction(scores, groups)
	perf <- performance(pred,"tpr", "fpr")
	
	idx <- .getSmallSeq(perf@x.values[[1]])
	plot( perf@x.values[[1]][idx], perf@y.values[[1]][idx], type="l", lwd=4, col="blue", xlab="FPR", ylab="TPR", xlim=xlim, main=main, cex.main=1.4, cex.lab=1.4, cex.axis=1.2 )
	for(i in 2:length(scoresList)) {
		idx <- .getSmallSeq(perf@x.values[[i]])
		lines( perf@x.values[[i]][idx], perf@y.values[[i]][idx], type="l", lwd=4, col=cols[i] )
	}
	legend("bottomright",legend=labels,lwd=4,col=cols,cex=1)
	invisible(list(tf=tf, scores=scores, groups=groups, keep=keep))
}



