
library(Repitools)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg18)
library(rtracklayer)

load("Rdata_created/orig_regs_gb.Rdata")
grps <- c("L=1 P=2","L=2 P=2","L=3 P=2","L=4 P=2","L=5 P=2")

# -----------------------------------------
# MBD-seq analyses
# -----------------------------------------

# ChIPDiff
cdpeaks <- read.table("peak_finders/ChIPDiff/LNCAP_vs_PREC.region")
cdGR <- GRanges(cdpeaks$V1, ranges=IRanges(cdpeaks$V2,cdpeaks$V3))
export(cdGR, "other/MBD_ChIPDiff.bed")
fo <- findOverlaps(cdGR, gb, select="first")
cdt <- table(regs[fo],cdpeaks$V4)

# RSEG
rspeaks <- read.table("peak_finders/rseg/P:L-domains.bed")
rsGR <- GRanges(rspeaks$V1, ranges=IRanges(rspeaks$V2,rspeaks$V3))
export(rsGR[grep("ENRICHED",rspeaks$V4)], "other/MBD_RSEG.bed")
fo <- findOverlaps(rsGR, gb, select="first")
rst <- table(regs[fo],rspeaks$V4)

load("Rdata_provided/DiffBind_analysis.Rdata")
fo <- findOverlaps(m.db.gr, gb, select="first")
dbt <- table(regs[fo],sign(values(m.db.gr)$Fold))
export(m.db.gr, "other/MBD_DiffBind.bed")
rm(m.db.gr)

load("Rdata_provided/DiffBind_analysis_nosubtract.Rdata")
fo <- findOverlaps(m.db.gr, gb, select="first")
dbt1 <- table(regs[fo],sign(values(m.db.gr)$Fold))

load("Rdata_created/MBD_abcdDNA_analysis.Rdata")
sig <- p.adjust( lrtadj$table$PValue, method="fdr") < .05
fo <- findOverlaps(gb.ABCD[sig], gb, select="first")
abt <- table(regs[fo],sign(lrtadj$table$logFC[sig]))
export(gb.ABCD[sig], "other/MBD_ABCD.bed")

# cdt - ChIPDiff
# rst - RSEG
# dbt - DiffBind
# abt - ABCD-DNA

pdf("Figures/Figure_4_RRPDs.pdf",h=5,w=15)
par(mfrow=c(1,3))
plot( (.5+cdt[,2])/cdt[,1], type="b", lwd=4, ylim=c(0.01,6), pch=19, log="y", xaxt="n", xlab="Relative CN state", ylab="Relative Peak Density", cex.lab=1.4, cex.axis=1.4)
lines((.5+rst[,3])/rst[,2],lwd=4,lty=1, col="orange", type="b", pch=19)
lines((.5+dbt[,2])/(dbt[,1]+.5),lwd=4,lty=1,col="grey", type="b", pch=19)
lines((.5+dbt1[,2])/(dbt1[,1]+.5),lwd=4,lty=1,col="grey30", type="b", pch=19)
lines(abt[,2]/abt[,1],lwd=4,lty=1,col="red", type="b", pch=19)
abline(h=1,col="grey",lty=3)
axis(1,labels=grps,at=1:length(grps),cex.axis=1.4)
legend("bottomright",c("ChIPDiff","RSEG","DiffBind-subtract","DiffBind-no subtract","ABCD-DNA"),lwd=4,lty=1,col=c("black","orange","grey","grey30","red"),cex=1.2)
title(main="MBDCap",cex.main=1.2)
mtext("A",3, cex=2, adj=-.1, padj=-.8)
#dev.off()



# -----------------------------------------
# H3K27me3-seq analyses
# -----------------------------------------

# ChIPDiff
cdpeaks <- read.table("peak_finders/ChIPDiff/LNCAP_vs_PREC_H3K27.region")
cdGR <- GRanges(cdpeaks$V1, ranges=IRanges(cdpeaks$V2,cdpeaks$V3))
fo <- findOverlaps(cdGR, gb, select="first")
cdt <- table(regs[fo],cdpeaks$V4)

# RSEG
rspeaks <- read.table("peak_finders/rseg/PrEC_H3K27me3:LNCaP_H3K27me3_subsetted-domains.bed")
rsGR <- GRanges(rspeaks$V1, ranges=IRanges(rspeaks$V2,rspeaks$V3))
fo <- findOverlaps(rsGR, gb, select="first")
rst <- table(regs[fo],rspeaks$V4)

load("Rdata_provided/DiffBind_analysis_K27.Rdata")
fo <- findOverlaps(m.db.gr, gb, select="first")
dbt <- table(regs[fo],sign(values(m.db.gr)$Fold))

load("Rdata_provided/DiffBind_analysis_K27_nosubtract.Rdata")
fo <- findOverlaps(m.db.gr, gb, select="first")
dbt1 <- table(regs[fo],sign(values(m.db.gr)$Fold))


load("Rdata_created/K27_abcdDNA_analysis.Rdata")
sig <- p.adjust( lrtadj$table$PValue, method="fdr") < .05
fo <- findOverlaps(gb.ABCD[sig], gb, select="first")
abt <- table(regs[fo],sign(lrtadj$table$logFC[sig]))

# cdt - ChIPDiff
# rst - RSEG
# dbt - DiffBind
# abt - ABCD-DNA

#pdf("Figures/Figure_4_RRPDs_H3K27.pdf",h=5,w=5)
plot( (.5+cdt[,2])/cdt[,1], type="b", lwd=4, ylim=c(0.01,15), log="y", xaxt="n", xlab="Relative CN state", ylab="Relative Peak Density", pch=19, cex.lab=1.4, cex.axis=1.4)
lines((.5+rst[,3])/rst[,2],lwd=4,lty=1, col="orange", type="b", pch=19)
lines((.5+dbt[,2])/(dbt[,1]+.5),lwd=4,lty=1,col="grey", type="b", pch=19)
lines((.5+dbt1[,2])/(dbt1[,1]+.5),lwd=4,lty=1,col="grey30", type="b", pch=19)
lines(abt[,2]/abt[,1],lwd=4,lty=1,col="red", type="b", pch=19)
abline(h=1,col="grey",lty=3)
axis(1,labels=grps,at=1:length(grps),cex.axis=1.4)
legend("bottomright",c("ChIPDiff","RSEG","DiffBind-subtract","DiffBind-no subtract","ABCD-DNA"),lwd=4,lty=1,col=c("black","orange","grey","grey30","red"),cex=1.2)
title(main="H3K27me3",cex.main=1.2)
mtext("B",3, cex=2, adj=-.1, padj=-.8)
#dev.off()


# -----------------------------------------
# H3K4me3-seq analyses
# -----------------------------------------

cdpeaks <- read.table("peak_finders/ChIPDiff/LNCAP_vs_PREC_H3K4.region")
cdGR <- GRanges(cdpeaks$V1, ranges=IRanges(cdpeaks$V2,cdpeaks$V3))
fo <- findOverlaps(cdGR, gb, select="first")
cdt <- table(regs[fo],cdpeaks$V4)

# RSEG
rspeaks <- read.table("peak_finders/rseg/PrEC_H3K4me3:LNCaP_H3K4me3-domains.bed")
rsGR <- GRanges(rspeaks$V1, ranges=IRanges(rspeaks$V2,rspeaks$V3))
fo <- findOverlaps(rsGR, gb, select="first")
rst <- table(regs[fo],rspeaks$V4)

load("Rdata_created/K4_abcdDNA_analysis.Rdata")
sig <- p.adjust( lrtadj$table$PValue, method="fdr") < .05
fo <- findOverlaps(gb.ABCD[sig], gb, select="first")
abt <- table(regs[fo],sign(lrtadj$table$logFC[sig]))

# cdt - ChIPDiff
# rst - RSEG
# dbt - DiffBind  # cannot do this one
# abt - ABCD-DNA

#pdf("Figures/Figure_4_RRPDs_H3K4.pdf",h=5,w=5)
plot( (.5+cdt[,2])/cdt[,1], type="b", lwd=4, ylim=c(0.001,9), log="y", xaxt="n", xlab="Relative CN state", ylab="Relative Peak Density", pch=19, cex.lab=1.4, cex.axis=1.4)
lines((.5+rst[,3])/rst[,2],lwd=4,lty=1, col="orange", type="b", pch=19)
lines(abt[,2]/abt[,1],lwd=4,lty=1,col="red", type="b", pch=19)
abline(h=1,col="grey",lty=3)
axis(1,labels=grps,at=1:length(grps),cex.axis=1.4)
legend("bottomright",c("ChIPDiff","RSEG","ABCD-DNA"),lwd=4,lty=1,col=c("black","orange","red"),cex=1.2)
title(main="H3K4me3",cex.main=1.2)
mtext("C",3, cex=2, adj=-.1, padj=-.8)
dev.off()


