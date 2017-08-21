
counts <- counts[k,]
gb <- gb[k]
regs <- regs[k]

grp <- regs
grp1 <- gsub(" P=2","",grp)

# calculate standard p-values / CN-aware p-values
# dataset-specific
d <- DGEList(counts=counts,group=sampGroup)
d <- calcNormFactors(d)

nm <- names(table(grp1))

# check on linearity
png(paste("Figures/Supplementary_Figure_3",letLabel,"_",samps,"_smear_by_CN.png",sep=""),w=length(nm)*400,h=400)
par(mfrow=c(1,length(nm)), mai=c(1.02,1.02,0.82,0.22))
f.by.cn <- rep(NA,length(nm))
s2 <- sum(d$counts[,ref])
s3 <- sum(d$counts[,sam])

n <- rep(NA,length(nm))


for(i in 1:length(nm)) {
	kk <- grp1==nm[i]
	map <- maPlot(d$counts[kk,ref]/s2, d$counts[kk,sam]/s3, normalize=FALSE, pch=19, cex=.3, ylim=c(-5,5), cex.main=3, cex.lab=3, cex.axis=3); grid();
	if(i==1) mtext(letLabel,3, cex=3, adj=-.1, padj=-.8)
	n[i] <- length(map$M)
	top <- max(100,round(.01*sum(!map$w)))
    #q <- quantile(map$A[!map$w], .95)
	o <- order(-map$A)[1:top]
	q <- max(min(map$A[o]),max(map$A[map$w]))
	cat(q,"\n")
	abline(v=q,col="blue"); 
        abline(h=median(map$M[map$A>q]),col="red",lwd=4); 
        title(paste("P=2 ",nm[i]," (",n[i],")",sep=""),cex.main=3)
	f.by.cn[i] <- median(map$M[map$A>q])
}
dev.off()

f <- exp(f.by.cn)
save(f,n, file=paste("Rdata_created/f.",samps,".Rdata",sep=""))


png(paste("other/NOTUSED_actual_normalization_",samps,".png",sep=""),w=1500,h=1000)
ref <- 1
nf <- rep(0,ncol(d))
par(mfrow=c(2,3))
# dataset-specific
for(i in 2:ncol(counts)) {
  cat(i,"\n")
  map <- maPlot(d$counts[grp1=="L=4",ref], d$counts[grp1=="L=4",i], normalize=TRUE, pch=19, cex=.3, ylim=c(-5,5)); grid();
  q <- quantile(map$A[!map$w],.99)
  nf[i] <- median( map$M[map$A > q] )
  abline(h=nf[i],col="red",lwd=4)
  abline(v=q,col="blue")
}
dev.off()

z <- exp(nf)
z <- z/(prod(z)^(1/length(z)))
d$samples$norm.factors <- z

# Fit the NB GLMs
design <- model.matrix(~d$samples$group)
fit <- glmFit(d, design, dispersion=0.05)
lrt <- glmLRT(fit, coef=2)
glmz <- -sign(lrt$table$logFC)*abs(qnorm(lrt$table$PValue/2))

# create matrix of CNs, then matrix of offsets
cnL <- as.numeric(gsub("L=","",gsub(" P=2","",regs)))
cn <- cbind(cnL/2,2) %*% m

o <- outer( rep(1,nrow(d)), getOffset(d)) + log(cn)

# Fit the CN-aware NB GLMs
fitadj <- glmFit(d, design, offset=o, dispersion=.05)
lrtadj <- glmLRT(fitadj, coef=2)
glmadjz <- -sign(lrtadj$table$logFC)*abs(qnorm(lrtadj$table$PValue/2))

gb.ABCD <- gb
save(fit,lrt,gb.ABCD,fitadj,lrtadj,file=paste("Rdata_created/",samps,"_abcdDNA_analysis.Rdata",sep=""))


