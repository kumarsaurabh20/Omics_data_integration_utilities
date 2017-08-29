### Testing the utility of an integrated analysis of copy number and transcriptomics datasets for inferring gene regulatory relationships
### Xin Yi Goh, Richard Newton*, Lorenz Wernisch, Rebecca Fitzgerald
### * richard.newton@mrc-bsu.cam.ac.uk

### R functions and example analysis

########## Example single dataset analysis ###########

if(FALSE){

	## Required 

	source("Code_S1.R")

	## Load example datasets - esophageal acgh and expression microarray data from the same samples
	## columns = samples, rows = probes, rownames = probenames.
	
	load("Data_S1.Rdata")

	## Function 'null.distb' generates the null distribution of partial correlations from a randomised dataset.

	nd  <- null.distb(eac.acgh, eac.expr)

	## Function 'find.pcors' calculates partial correlations between acgh and expression profiles

	ip <- find.pcors(eac.acgh, eac.expr)

	## Function 'find.g1s' generates list of probes ranked by significance of partial 
	## correlation between own acgh and expression profiles
	## ip = output from 'find.pcors', nd = output from 'null.distb'

	g1s <- find.g1s(ip, nd)

	## Function 'find.g2s' generates list of potential target probes for a potential regulatory gene (g1)
	## g1 = name of potential regulatory gene, ip = output from 'find.pcors', 
	## nd = output from 'null.distb', alt = 1 or -1 for positive or negative correlation.

	g2s.erbb2.pos <- find.g2s(g1="ERBB2", ip, nd, 1)
	g2s.erbb2.neg <- find.g2s(g1="ERBB2", ip, nd, -1)

	g2s.fgfr2.pos <- find.g2s(g1="FGFR2", ip, nd, 1)
	g2s.fgfr2.neg <- find.g2s(g1="FGFR2", ip, nd, -1)

}


##### Single dataset analysis functions #####

library(GeneNet)

null.distb <- function(acgh, expr){
    N <- nrow(acgh)
    comb<-rbind(expr[sample(1:N),], acgh[sample(1:N),])  
    inf.pcor <- ggm.estimate.pcor(t(comb))                    
    rm(comb); gc()
    ip <- inf.pcor[1:N,(N+1):(N*2)]
    ip <- sample(ip, 1e06)
    rm(inf.pcor); gc()
    nd <- c(mean(ip), sd(ip))
}

find.pcors <- function(acgh, expr){
    N <- nrow(acgh)
    comb <- rbind(expr, acgh)                         
    nmes <- rownames(comb)
    inf.pcor <- ggm.estimate.pcor(t(comb))            
    rm(comb); gc()
    ip <- inf.pcor[1:N,(N+1):(N*2)]
    rm(inf.pcor); gc()
    rownames(ip)<-rownames(expr)
    colnames(ip) <- rownames(expr)
    ip
}

find.g1s <- function(ip, nd){  
    ip.g1 <- diag(ip)                           
    g1s <- sapply(ip.g1, function(x) pnorm(x, nd[1], nd[2], lower.tail=F))  ## p-value for pcors from null distribution
    names(g1s) <- rownames(ip)
    g1s <- g1s[order(g1s)]
    g1s <- p.adjust(g1s, method="fdr")                                          
}

find.g2s <- function(g1, ip, nd, alt){
    g1.ind <- which(colnames(ip)==g1)
    g1c.g2e <- sapply(ip[,g1.ind], function(x) pnorm(sign(alt)*x, nd[1], nd[2], lower.tail=F))
    g2c.g1e <- sapply(ip[g1.ind,], function(x) pnorm(sign(alt)*x, nd[1], nd[2], lower.tail=F))
    g1c.g2e <- p.adjust(g1c.g2e, method="fdr"); g2c.g1e <- p.adjust(g2c.g1e, method="fdr")
    g2s <- data.frame(g2=rownames(ip), g1c.g2e, g2c.g1e)
    g2s <- g2s[order(g2s$g1c.g2e),]
}


##### Multiple dataset analysis  functions ######

## m.acgh and m.expr = lists of matched acgh/expression datasets

library(survcomp)

gs.in.common <- function(m.acgh){
    gs.cmn <- rownames(m.acgh[[1]])
    for(i in 2:length(m.acgh)){
        gs.cmn <- intersect(gs.cmn, rownames(m.acgh[[i]]))
    }
    gs.cmn
}

## Pearson correlations ##

null.distb.m <- function(m.acgh, m.expr, its, verbose){
    M <- length(m.acgh)
    corr.pv <- array(NA, c(its,M))
    for(j in 1:M){
        if(verbose){print(j)}
        acgh <- m.acgh[[j]]; expr <- m.expr[[j]]; N <- nrow(acgh)
        for(i in 1:its){
            ind.a <- sample(1:N,1)
            ind.e <- sample(1:N,1)
            ac <- acgh[ind.a,]; ex <- expr[ind.e,]
            if(length(which(!is.na(ex*ac)))>2){
                if(sd(ac, na.rm=T)>0 && sd(ex, na.rm=T)>0){
                    pv <- cor.test(ac, ex, alternative="g")$p.value
                    if(pv == 0){
                        pv <- 2e-16
                    }
                    corr.pv[i,j] <- pv
                }
            }
        }
    }
    comb.pv.rand <- c(0, apply(corr.pv, 1, combine.test, na.rm=T))
}
        
find.g1s.m <- function(gs.cmn, m.acgh, m.expr, comb.pv.rand){
    M <- length(m.acgh)
    N <- length(gs.cmn)
    corr.pv <- array(NA, c(N,M))
    inds.na <- NULL
    for(j in 1:M){
        acgh <- m.acgh[[j]]; expr <- m.expr[[j]]
        for(i in 1:N){
            ind <- which(rownames(acgh) == gs.cmn[i])
            ac <- acgh[ind,]; ex <- expr[ind,]
            if(length(which(!is.na(ex*ac)))>2){
                pv <- cor.test(ac, ex, alternative="g")$p.value
                if(pv == 0){
                    pv <- 2e-16
                }
                corr.pv[i,j] <- pv
            }else{
                inds.na <- c(inds.na, i)
            }
        }
    }
    rownames(corr.pv) <- gs.cmn
    inds.na <- unique(inds.na); corr.pv <- corr.pv[-inds.na,]
    comb.pv <- apply(corr.pv, 1, combine.test, na.rm=T)
    perm.pv <- sapply(comb.pv, function(x) length(which(comb.pv.rand<=x))/length(comb.pv.rand))
    perm.pv.adj <- p.adjust(perm.pv, method="fdr")
    sams <- corr.pv; sams[which(sams<=0.05)] <- "X"; sams[which(sams != "X")] <- ""
    num <- apply(sams, 1, function(x) length(which(x=="X")))
    g1s.m <- data.frame(g1=rownames(corr.pv), pv.adj=perm.pv.adj, num=num, sams)
    g1s.m <- g1s.m[order(g1s.m[,3], 1-g1s.m[,2], decreasing=T),]
}

find.g2s.m <- function(g1, m.acgh, m.expr, g1s.m, alt, comb.pv.rand, verbose=F){
    M <- length(m.acgh)
    g2s.unq <- unique(unlist(lapply(m.acgh, rownames))) 
    N <- length(g2s.unq)
    ac1ex2pv<- matrix(NA, N, M)
    ac2ex1pv<- matrix(NA, N, M)
    for(j in 1:M){
        if(verbose){print(j)}
        acgh <- m.acgh[[j]]; expr <- m.expr[[j]]
        g1.ind <- which(rownames(acgh)==g1)
        ac.g1.s <- acgh[g1.ind,]; ex.g1.s <- expr[g1.ind,]
        inds.g1 <- intersect(which(!is.na(ac.g1.s)), which(!is.na(ex.g1.s)))
        for(i in 1:N){
            g2.ind <- which(rownames(acgh) == g2s.unq[i])
            if(length(g2.ind)>0){
                ac.g2 <- acgh[g2.ind,]; ex.g2 <- expr[g2.ind,] 
                inds.g2 <- intersect(which(!is.na(ac.g2)), which(!is.na(ex.g2)))
                inds <- intersect(inds.g1, inds.g2)
                if(length(inds)>2){
                    ac.g1 <- ac.g1.s[inds]; ex.g1 <- ex.g1.s[inds]
                    ac.g2 <- ac.g2[inds];   ex.g2 <- ex.g2[inds]
                    if(sd(ac.g1)>0 && sd(ex.g1) > 0 && sd(ac.g2) >0 && sd(ex.g2)>0){
                        ac1ex2pv[i,j] <-  cor.test(ac.g1, ex.g2, alternative=alt)$p.value
                        ac2ex1pv[i,j] <- cor.test(ac.g2, ex.g1,  alternative=alt)$p.value
                    }
                }
            }
        }
    }
    if(verbose){print("Calculating p-values ...")}
    g1.amp.inds <- which(g1s.m[which(g1s.m$g1 == g1), 4:(M+3)] == "X")  
    comb.pv <- apply(ac1ex2pv[,g1.amp.inds], 1, combine.test, na.rm=T)
    perm.pv <- sapply(comb.pv, function(x) length(which(comb.pv.rand<=x))/length(comb.pv.rand))
    perm.pv.adj <- p.adjust(perm.pv, method="fdr")
    sams <- ac1ex2pv; sams[which(sams<=0.05)] <- "X"
    sams[which(sams != "X")] <- ""; sams[,setdiff(1:M,g1.amp.inds)] <- ""; sams[which(is.na(sams))]<-""
    num <- apply(sams, 1, function(x) length(which(x=="X")))
    comb.pv.r <- apply(ac2ex1pv[,g1.amp.inds], 1, combine.test, na.rm=T)
    perm.pv.r <- sapply(comb.pv.r, function(x) length(which(comb.pv.rand<=x))/length(comb.pv.rand))
    perm.pv.r.adj <- p.adjust(perm.pv.r, method="fdr")
    g2s.m <- data.frame(g2=g2s.unq, pv.adj=perm.pv.adj, pv.adj.rev=perm.pv.r.adj, num=num, sams)
    g2s.m <- g2s.m[order(g2s.m[,4], 1-g2s.m[,2], decreasing=T),]
}



## Partial correlations ##

library(GeneNet)

null.distb.m.pc.rnd <- function(m.acgh, m.expr, len=1e06, verbose=T){
    M <- length(m.acgh)
    corr.pv <- array(NA, c(len,M))
    nd.pc.rnd <-array(NA, c(2,M))
    for(j in 1:M){
        if(verbose){print(j)}
        acgh <- m.acgh[[j]]; expr <- m.expr[[j]]; N <- nrow(acgh)
        N <- nrow(acgh)
        C <- ncol(acgh)
        comb<-rbind(expr[sample(1:N),sample(1:C)], acgh[sample(1:N),sample(1:C)])  
        inf.pcor <- ggm.estimate.pcor(t(comb))                    
        rm(comb); gc()
        ip <- inf.pcor[1:N,(N+1):(N*2)]
        rm(inf.pcor); gc()
        ip <- as.vector(sample(ip, len))
        nd.pc.rnd[1,j] <- mean(ip)
        nd.pc.rnd[2,j] <- sd(ip)
    }
    nd.pc.rnd
}

null.distb.m.pc <- function(m.acgh, m.expr, nd.pc.rnd, len=1e06, verbose=T){
    M <- length(m.acgh)
    corr.pv <- array(NA, c(len,M))
    for(j in 1:M){
        if(verbose){print(j)}
        acgh <- m.acgh[[j]]; expr <- m.expr[[j]]; N <- nrow(acgh)
        N <- nrow(acgh)
        comb<-rbind(expr[sample(1:N),], acgh[sample(1:N),])  
        inf.pcor <- ggm.estimate.pcor(t(comb))                    
        rm(comb); gc()
        ip <- inf.pcor[1:N,(N+1):(N*2)]
        rm(inf.pcor); gc()
        ip <- as.vector(sample(ip, len))
        corr.pv[,j] <- sapply(ip, function(x) pnorm(x, nd.pc.rnd[1,j], nd.pc.rnd[2,j], lower.tail=F))
    }
    comb.pv.rand <- c(0, apply(corr.pv, 1, combine.test, na.rm=T))
}

        
find.g1s.m.pc <- function(gs.cmn, m.acgh, m.expr, comb.pv.rand, nd.pc.rnd, verbose=T){
    M <- length(m.acgh)
    G <- length(gs.cmn)
    corr.pv <- array(NA, c(G,M))
    for(j in 1:M){
        acgh <- m.acgh[[j]]; expr <- m.expr[[j]]
        if(verbose){print(j)}
        N <- nrow(acgh)
        comb<-rbind(expr,acgh)  
        inf.pcor <- ggm.estimate.pcor(t(comb))                    
        rm(comb); gc()
        ip <- diag(inf.pcor[1:N,(N+1):(N*2)])
        rm(inf.pcor); gc()
        tt <- sapply(ip, function(x) pnorm(x, nd.pc.rnd[1,j], nd.pc.rnd[2,j], lower.tail=F))
        names(tt) <- rownames(acgh)
        inds <- which(names(tt) %in% gs.cmn)
        tt <-tt[inds]
        tt <- tt[order(names(tt))]
        corr.pv[,j] <- tt
        rownames(corr.pv) <- names(tt)
    }
    comb.pv <- apply(corr.pv, 1, combine.test, na.rm=T)
    perm.pv <- sapply(comb.pv, function(x) length(which(comb.pv.rand<=x))/length(comb.pv.rand))
    perm.pv.adj <- p.adjust(perm.pv, method="fdr")
    sams <- corr.pv; sams[which(sams<=0.05)] <- "X"; sams[which(sams != "X")] <- ""
    num <- apply(sams, 1, function(x) length(which(x=="X")))
    g1s.m <- data.frame(g1=rownames(corr.pv), pv.adj=perm.pv.adj, num=num, sams)
    g1s.m <- g1s.m[order(g1s.m[,3], 1-g1s.m[,2], decreasing=T),]
}

find.g2s.m.pc <- function(g1, m.acgh, m.expr, g1s.amp, alt, comb.pv.rand, nd.pc.rnd, verbose=F){
    M <- length(m.acgh)
    g2s.unq <- unique(unlist(lapply(m.acgh, rownames)))
    g2s.unq <- g2s.unq[order(g2s.unq)]
    G <- length(g2s.unq)
    ac1ex2pv<- matrix(NA, G, M)
    rownames(ac1ex2pv) <- g2s.unq
    for(j in 1:M){
        if(verbose){print(j)}
        acgh <- m.acgh[[j]]; expr <- m.expr[[j]]
        if(g1 %in% rownames(acgh)){
            N <- nrow(acgh)
            comb<-rbind(expr,acgh)  
            inf.pcor <- ggm.estimate.pcor(t(comb))                   
            rm(comb); gc()
            ip <- inf.pcor[1:N,(N+1):(N*2)]
            rm(inf.pcor); gc()
            g1.ind <- which(colnames(ip)==g1)
            g1c.g2e <- sapply(ip[,g1.ind], function(x) pnorm(sign(alt)*x, nd.pc.rnd[1,j], nd.pc.rnd[2,j], lower.tail=F))
            names(g1c.g2e) <- rownames(ip)
            g1c.g2e <- g1c.g2e[order(names(g1c.g2e))]
            inds <- which(g2s.unq %in% names(g1c.g2e))
            ac1ex2pv[inds, j] <- g1c.g2e
        }
    }
    g1.amp.inds <- g1s.amp 
    ac1ex2pv <- ac1ex2pv[,g1.amp.inds]
    comb.pv <- apply(ac1ex2pv, 1, combine.test, na.rm=T)
    perm.pv <- sapply(comb.pv, function(x) length(which(comb.pv.rand<=x))/length(comb.pv.rand))
    perm.pv.adj <- p.adjust(perm.pv, method="fdr")
    names(perm.pv.adj) <- g2s.unq
    perm.pv.adj <- perm.pv.adj[order(perm.pv.adj)]
}

