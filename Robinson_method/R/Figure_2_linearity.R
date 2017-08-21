
library(Repitools)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg18)
library(edgeR)

# read MBD data, count, normalize, do ABCD-DNA
load("Rdata_created/orig_regs_gb.Rdata")
load("Rdata_created/counts.Rdata")
rs <- rowSums(counts[,5:6])  # read depth in SssI
k <- rs>10
counts <- counts[,1:4]
sampGroup <- rep(c("LNCaP","PrEC"),each=2)
samps <- "MBD"
ref <- 3
sam <- 1
m <- matrix(rep(c(1,0,1),c(2,4,2)),nrow=2,byrow=TRUE)
letLabel <- "A"
source("R/run_ABCD-DNA.R", verbose=TRUE)


# read H3K27me3 data, count, normalize, do ABCD-DNA
load("Rdata_created/orig_regs_gb.Rdata")
load("Rdata_provided/H3K27me3_grl.Rdata")
counts <- annotationBlocksCounts(grl, gb, seq.len=300)
k <- rowSums(counts[,1:4]) > 15 | rowSums(counts[,5:7]) > 15
sampGroup <- rep(c("LNCaP","PrEC"),c(4,3))
samps <- "K27"
ref <- 5
sam <- 1
m <- matrix(rep(c(1,0,1),c(4,7,3)),nrow=2,byrow=TRUE)
letLabel <- "B"
source("R/run_ABCD-DNA.R", verbose=TRUE)


# read H3K4me3 data, count, normalize, do ABCD-DNA
load("Rdata_created/orig_regs_gb.Rdata")
load("Rdata_provided/H3K4me3_grl.Rdata")
counts <- annotationBlocksCounts(grl, gb, seq.len=300)
k <- rowSums(counts) > 10
sampGroup <- c("L","P")
samps <- "K4"
ref <- 2
sam <- 1
m <- matrix(rep(c(1,0,1),c(1,2,1)),nrow=2,byrow=TRUE)
letLabel <- "C"
source("R/run_ABCD-DNA.R", verbose=TRUE)


load("Rdata_created/f.K27.Rdata")
f.k27 <- f
load("Rdata_created/f.K4.Rdata")
f.k4 <- f
load("Rdata_created/f.MBD.Rdata")
f.mbd <- f

fs <- cbind(f.mbd, f.k27, f.k4)
cols <- c("red","blue","black")

pdf("Figures/Figure_2_linearity.pdf",h=5,w=5)
matplot( t(t(fs) / fs[4,]*4), lwd=4, type="b", lty=c(1,1,1) ,pch=20, col=cols, xlab="LNCaP copy number (PrEC copy=2)", ylab="Median fold-change (scaled to 4)", cex.lab=1.3)
#abline(0,1,col="grey",lty=3)
grid(col="grey")
legend("bottomright",c("MBDCap-seq","H3K27me3-seq","H3K4me3-seq"),lty=1,lwd=3,col=cols)
dev.off()
