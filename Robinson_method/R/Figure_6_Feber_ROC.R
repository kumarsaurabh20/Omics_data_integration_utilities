
library(limma)
library(Repitools)
library(GenomicRanges)
library(edgeR)
library(BSgenome.Hsapiens.UCSC.hg18)


load("Rdata_provided/m27kGR_wextras.Rdata")
load("Rdata_provided/beck.counts.Rdata")
load("Rdata_provided/batavgs.Rdata")
source("R/functions.R")

v <- values(m27kGR)

d <- DGEList(counts=counts, group=colnames(counts))
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)


# most prominent state
k <- v$cnN==3 & v$cnM==3 & v$cnS==3
k <- k & !is.na(k)

png("other/NOTUSED_actual_normalization_feber.png",w=1000,h=500)
ref <- 1
nf <- rep(0,ncol(d))
cs <- colSums(d$counts)
par(mfrow=c(1,2))
for(i in setdiff(1:ncol(d),ref)) {
	map <- maPlot(d$counts[k,ref]/cs[ref], d$counts[k,i]/cs[i], normalize=TRUE, pch=19, cex=.3, ylim=c(-5,5)); grid();
	q <- quantile(map$A[!map$w],.95)
	nf[i] <- median( map$M[map$A > q] )
	abline(h=nf[i],col="red",lwd=4)
	abline(v=q,col="blue")
}
dev.off()

z <- exp(-nf)
z <- z/(prod(z)^(1/length(z)))
d$samples$norm.factors <- z


newCh <- paste("chr",as.character(seqnames(m27kGR)),sep="")
m27kE <- GRanges(seqnames=newCh,ranges(m27kGR))
m27kE <- resize(m27kE,600,fix="center")
#cpgdens <- cpgDensityCalc(m27kE, organism=Hsapiens, w.function="none")




# 1. M-S (malignant-schwann = cancer-normal) ----------
design <- model.matrix(~d$samples$group)
fit <- glmFit(d, design, offset=getOffset(d), dispersion=.05)
lrt <- glmLRT(fit, coef=3)
glmz <- -sign(lrt$table$logFC)*abs(qnorm(lrt$table$PValue/2))

o <- outer( rep(1,nrow(d)), getOffset(d) ) + log( cbind(v$cnM,v$cnN,v$cnS) )
ww <- rowSums(is.na(o)) == 0
fitadj <- glmFit(d[ww,], design, offset=o[ww,], dispersion=.05)
lrtadj <- glmLRT(fitadj, coef=3)
glmadjz <- rep(NA, length(glmz))
glmadjz[ww] <- -sign(lrtadj$table$logFC)*abs(qnorm(lrtadj$table$PValue/2))


scoresList <- list(naivez=glmz,batman=battab[,2]-battab[,3],cnadjz=glmadjz)

#fn <- "Figures/Figure_6a_ROC_cancer-normal.pdf"
fn <- "Figures/Figure_6_ROC_Feber.pdf"
pdf(fn,h=4,w=12)
par(mfrow=c(1,3))
makeROCCurve.overall(md=v$Mb-v$Sb,main="cancer-normal",symmetric=FALSE)
mtext("A",3, cex=2, adj=-.1, padj=-.8)
#dev.off()



# 2. M-N (malignant-benign = cancer-benign) ----------
fit <- glmFit(d, design, offset=getOffset(d), dispersion=.05)
lrt <- glmLRT(fit, coef=2)
glmz <- -sign(lrt$table$logFC)*abs(qnorm(lrt$table$PValue/2))

o <- outer( rep(1,nrow(d)), getOffset(d) ) + log( cbind(v$cnM,v$cnN,v$cnS) )
ww <- rowSums(is.na(o)) == 0
fitadj <- glmFit(d[ww,], design, offset=o[ww,], dispersion=.05)
lrtadj <- glmLRT(fitadj, coef=2)
glmadjz <- rep(NA, length(glmz))
glmadjz[ww] <- -sign(lrtadj$table$logFC)*abs(qnorm(lrtadj$table$PValue/2))


scoresList <- list(naivez=glmz,batman=battab[,2]-battab[,1],cnadjz=glmadjz)

#fn <- "Figures/Figure_6b_ROC_cancer-benign.pdf"
#pdf(fn,h=5,w=5)
makeROCCurve.overall(md=v$Mb-v$Nb,main="cancer-benign",symmetric=FALSE)
mtext("B",3, cex=2, adj=-.1, padj=-.8)
#dev.off()




# 3. N-S (neurofibroma-schwann = benign-normal) ----------
grp <- gsub("GSM[0-9][0-9][0-9][0-9][0-9][0-9]_","",d$samples$group)
design <- model.matrix(~grp)
fit <- glmFit(d, design, offset=getOffset(d), dispersion=.05)
lrt <- glmLRT(fit, coef=3)
glmz <- -sign(lrt$table$logFC)*abs(qnorm(lrt$table$PValue/2))

o <- outer( rep(1,nrow(d)), getOffset(d) ) + log( cbind(v$cnM,v$cnN,v$cnS) )
ww <- rowSums(is.na(o)) == 0
fitadj <- glmFit(d[ww,], design, offset=o[ww,], dispersion=.05)
lrtadj <- glmLRT(fitadj, contrast=c(0,1,-1))
glmadjz <- rep(NA, length(glmz))
glmadjz[ww] <- -sign(lrtadj$table$logFC)*abs(qnorm(lrtadj$table$PValue/2))


scoresList <- list(naivez=glmz,batman=battab[,1]-battab[,3],cnadjz=glmadjz)

#fn <- "Figures/Figure_6c_ROC_benign-normal.pdf"
#pdf(fn,h=5,w=5)
makeROCCurve.overall(md=v$Nb-v$Sb,main="benign-normal",symmetric=FALSE)
mtext("C",3, cex=2, adj=-.1, padj=-.8)
dev.off()


