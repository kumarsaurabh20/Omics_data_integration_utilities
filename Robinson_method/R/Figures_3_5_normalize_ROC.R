
library(Repitools)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg18)
library(edgeR)

load("Rdata_provided/illumina_450k_lncap_prec.Rdata")
load("Rdata_created/counts.Rdata")
load("Rdata_created/orig_regs_gb.Rdata")
source("R/functions.R")

# widdle it down a bit
rs <- rowSums(counts[,5:6])  # read depth in SssI
k <- rs>10

counts <- counts[k,]
gb <- gb[k,]
regs <- regs[k]

grp <- regs
grp1 <- gsub(" P=2","",grp)

# calculate standard p-values / CN-aware p-values
d <- DGEList(counts=counts[,1:4],group=rep(c("LNCaP","PrEC"),each=2))
d <- calcNormFactors(d)

png("other/NOTUSED_actual_normalization.png",w=1500,h=500)
ref <- 1
nf <- rep(0,ncol(d))
par(mfrow=c(1,3))
for(i in 2:4) {
  map <- maPlot(d$counts[grp1=="L=4",ref], d$counts[grp1=="L=4",i], normalize=TRUE, pch=19, cex=.3, ylim=c(-5,5), xlim=c(-22.5,-16)); grid();
  q <- quantile(map$A[!map$w],.99)
  nf[i] <- median( map$M[map$A > q] )
  abline(h=nf[i],col="red",lwd=4)
  abline(v=q,col="blue")
}
dev.off()

z <- exp(nf)
z <- z/(prod(z)^(1/length(z)))
d$samples$norm.factors <- z


ref <- 3
#
#png("Figures/Figure_3_normalize_figure.png",w=2400,h=1400)
pdf("Figures/Figure_3_normalize_figure.pdf",w=8,h=4)
par(mfrow=c(1,2),mai=c(.8,.8,.8,.2))
mains <- c("LNCaP vs. PrEC",NA,NA,"PrEC vs. PrEC")
labels <- c("B",NA,NA,"A")
for(i in c(4,1)) {
	map <- maPlot(d$counts[grp1=="L=4",ref], d$counts[grp1=="L=4",i],
                  normalize=TRUE, plot.it=FALSE)
	q <- quantile(map$A[!map$w],.99)
	m <- median( map$M[map$A > q] )
    
    s <- sample(length(map$w),2e4)
    plot(map$A[s], map$M[s], col=c("black","orange")[map$w[s]+1],
         pch=19, cex=.3, ylim=c(-5,5), xlim=c(-24,-15), 
         cex.lab=1, cex.axis=1, xlab="A", ylab="M"); grid();

	abline(h=m,col="red",lwd=4)
	abline(v=q,col="blue",lty=1, lwd=2)
	title(main=mains[i],cex.main=1.2)
	mtext(labels[i],3, cex=1.4, adj=-.1, padj=-.8)

}
dev.off()




cpgdens <- cpgDensityCalc(gb, organism=Hsapiens, w.function="none")
cpgcut <- cut(cpgdens[grp1=="L=4"], quantile(cpgdens,p=(0:3)/3), include.lowest=TRUE)

png("other/NOTUSED_plotsmear")
ps <- plotSmear( d[grp1=="L=4",], pair=c("PrEC","LNCaP"), ylim=c(-5,5), xlim=c(-24,-16) )
dev.off()


library(lattice)
png("Figures/Supplementary_Figure_4_normalize_by_cpgdens.png",h=600,w=1800)
xyplot( ps$M~ps$A | cpgcut, ylim=c(-5.5,5.5), xlim=c(-24.5,-15), panel=panelFunction, xlab="A", ylab="M")
dev.off()


# Fit the NB GLMs
design <- model.matrix(~d$samples$group)
d <- estimateGLMCommonDisp(d, design)
fit <- glmFit(d, design, dispersion=d$common.dispersion)
lrt <- glmLRT(fit, coef=2)
glmz <- -sign(lrt$table$logFC)*abs(qnorm(lrt$table$PValue/2))


# create matrix of CNs, then matrix of offsets
cnL <- as.numeric(gsub("L=","",gsub(" P=2","",regs)))
m <- matrix(rep(c(1,0,1),c(2,4,2)),nrow=2,byrow=TRUE)
cn <- cbind(cnL/2,2) %*% m

o <- outer( rep(1,nrow(d)), getOffset(d)) + log(cn)

# Fit the CN-aware NB GLMs
fitadj <- glmFit(d, design, offset=o, dispersion=d$common.dispersion)
lrtadj <- glmLRT(fitadj, coef=2)
glmadjz1 <- glmadjz <- -sign(lrtadj$table$logFC)*abs(qnorm(lrtadj$table$PValue/2))


# Fit CN-aware NB GLMs, but use segmented gDNA-seq instead
load("Rdata_provided/gDNAseq_CN_LNCAP_PREC_30kb.Rdata")

foL <- findOverlaps(gb, acn3@adj.CN.seg[[2]], select="first")
foP <- findOverlaps(gb, acn3@adj.CN.seg[[1]], select="first")

cnL <- values(acn3@adj.CN.seg[[2]])$CN[foL]
cnP <- values(acn3@adj.CN.seg[[1]])$CN[foP]

md <- median( cnL/cnP, na.rm=TRUE )
cn <- cbind(cnL/md,cnP) %*% m

o <- outer( rep(1,nrow(d)), getOffset(d)) + log(cn)

k <- rowSums( is.na(o) ) == 0

fitadj1 <- glmFit(d[k,], design, offset=o[k,], dispersion=d$common.dispersion)
lrtadj1 <- glmLRT(fitadj1, coef=2)
glmadjz1[k] <- -sign(lrtadj1$table$logFC)*abs(qnorm(lrtadj1$table$PValue/2))
glmadjz1[!k] <- NA


load("Rdata_provided/DiffBind_analysis.Rdata")
fo <- findOverlaps(gb, m.db.all.gr, select="first")
v <- values(m.db.all.gr)
dbz <- sign(v$Fold[fo])*abs(qnorm(v$p.value[fo]/2))

load("Rdata_provided/DiffBind_analysis_nosubtract.Rdata")
fo <- findOverlaps(gb, m.db.all.gr, select="first")
v <- values(m.db.all.gr)
dbz1 <- sign(v$Fold[fo])*abs(qnorm(v$p.value[fo]/2))


# add rseg in the mix
root <- "peak_finders/rseg/roc/"
f <- dir(root,"domains",recursive=TRUE,full=TRUE)[1:40]
names(f) <- gsub(root,"",f)
tabs <- lapply(f, read.table)
grs <- lapply(tabs, function(u) GRanges(u$V1,IRanges(u$V2,u$V3),grp=u$V4))
grs <- GRangesList(grs)

scores <- as.numeric(gsub("/","",gsub("/L_sorted:P_sorted-domains.bed","",names(f))))

z <- matrix(NA, nrow=length(gb), ncol=length(scores))
colnames(z) <- scores
for(i in 1:ncol(z)) {
	cat(names(grs)[i],"\n")
	gr <- grs[[i]]
	fo <- findOverlaps(gb, gr, select="first")
	z[,i] <- as.character(values(gr)$grp[fo])
}

ze <- matrix(z %in% c("SAMPLE-I-ENRICHED","SAMPLE-II-ENRICHED"),ncol=length(scores))

ind <- 1:length(scores)
sx <- scores[apply(ze,1,function(u) max(ind[u]))]
sx[is.na(sx)] <- 0

library(VennDiagram)

cutoff <- .01
w1 <- p.adjust(lrt$table$PValue, method="fdr") < cutoff
w2 <- p.adjust(lrtadj$table$PValue, method="fdr") < cutoff

pdf("Figures/Supplementary_Figure_6_VennDiagrams.pdf",w=9,h=5)
z <- venn.diagram(list(Naive=which(w1),"CNV-Aware"=which(w2)), main="Overall", filename=NULL, fill=c("red","blue"), cex=2, main.cex=2, cat.cex=2);
grid.draw(z)
for(i in 1:5) {
	grid.newpage()
	k <- grp1==paste("L=",i,sep="")
	z <- venn.diagram(list(Naive=which(w1[k]),"CNV-Aware"=which(w2[k])), main=grp1[k][1], filename=NULL, fill=c("red","blue"), cex=2, main.cex=2, cat.cex=2);
	grid.draw(z)
}
dev.off()


cngrp <- names(table(grp))[-1]
#cngrp <- names(table(grp))

# marry up the 450k array data to the MBD-seq
## --------------------------------------
load("Rdata_provided/illumina_450k_lncap_prec.Rdata")
fo <- findOverlaps(mGR, gb)
#inds <- split( fo@matchMatrix[,1], fo@matchMatrix[,2] )  # old BioC
inds <- split( fo@queryHits, fo@subjectHits )  # new BioC
v <- as.data.frame(values(mGR))
betaavg <- sapply(inds, function(u) { cat("."); colMeans(v[u,,drop=FALSE])} )
betaavg <- t(betaavg)
betavals <- matrix(NA, nrow=nrow(counts),ncol=2)
ind <- as.integer(names(inds))
betavals[ind,] <- betaavg


# partition 450k regions into true/false
high.diff <- .3
low.diff <- .1

tf <- rep(NA, nrow(betavals))
md <- betavals[,2]-betavals[,1]
tf[abs(md) > high.diff] <- 1
tf[abs(md) < low.diff] <- 0
tf[is.na(tf)] <- -99
keep <- which(tf >= 0 & !is.na(tf))

table(tf,sign(md),regs)
sho <- tf>=0 & grp %in% cngrp

md[tf==0] <- 0

pdf("Figures/Supplementary_Figure_7_Zscores_by_CNV_truth.pdf",w=11,h=8.5)
par(mfrow=c(2,1))
boxplot( glmz ~ sign(md)+grp1, subset=sho, cex.axis=.8, las=2, outline=FALSE, ylim=c(-6,6), col=c("lightgreen","lightblue","lightgreen"), main="Naive Z-scores"); abline(h=0)
legend("topleft",c("Non-Differential","Differential"),pch=15,col=c("lightblue","lightgreen"))
abline(v=(1:3)*3+.5,col="orange")
boxplot( glmadjz ~ sign(md)+grp1, subset=sho, cex.axis=.8, las=2, outline=FALSE, ylim=c(-6,6), col=c("lightgreen","lightblue","lightgreen"), main="CNV-aware Z-scores"); abline(h=0)
abline(v=(1:3)*3+.5,col="orange")
dev.off()


rm(tf,md)


# overall ROC curves
library(ROCR)

scoresList <- list(naive=glmz,adjusted=glmadjz,gdna=glmadjz1,rseg=sx,db=dbz, db1=dbz1)

rocCols=c("blue","red","salmon","orange","grey","grey30")
rocLabels=c("Naive","ABCD-DNA (SNP 6.0)","ABCD-DNA (gDNA-seq)","RSEG","DiffBind-subtract","DiffBind-no subtract")

pdf("Figures/Figure_5_ROC_lncap_prec.pdf",h=8,w=8)
par(mfrow=c(2,2),mai=c(.8,.7,.6,.1)); s <- makeROCCurve(high.diff=.4,low.diff=.1,xlim=c(0,1), symmetric=TRUE, labels=rocLabels, cols=rocCols, min.cpgw=7)
dev.off()

pdf("Figures/Supplementary_Figure_8_ROC_nonsymmetric.pdf",h=8,w=8)
par(mfrow=c(2,2)); s <- makeROCCurve(high.diff=.4,low.diff=.1,xlim=c(0,1), symmetric=TRUE, labels=rocLabels, cols=rocCols, min.cpgw=7)
dev.off()

#pdf("Figures/Supplementary_Figure__8_roc_symmetric_strict.pdf",h=8,w=8)
#par(mfrow=c(2,2)); s <- makeROCCurve(high.diff=.5,low.diff=.05,xlim=c(0,1), symmetric=TRUE, labels=rocLabels, cols=rocCols, min.cpgw=7)
#dev.off()

