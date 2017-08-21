## pre-processing of PICNIC SNP array data
## --------------------------------------
library(Repitools)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg18)

# read in CN info
# partition genome into copy number regions
# -----------------------------------------
cnL <- read.csv("PICNIC/segments__feature.TXT_LNCaP.csv", header=FALSE)
cnP <- read.csv("PICNIC/segments__feature.TXT_PrEC.csv", header=FALSE)

picnic2GR <- function(x) {
	ch <- paste("chr",gsub("24","Y",gsub("23","X",x$V4)),sep="")
	GRanges( seqnames=ch,ranges=IRanges(start=x$V5,end=x$V6), picnic=DataFrame(cn=x$V8))
}
cnPgr <- picnic2GR(cnP)
cnLgr <- picnic2GR(cnL)

# make genome-wide bins and attribute copy number state to each
gb <- genomeBlocks(seqlengths(Hsapiens)[1:24], width=500)
foP <- findOverlaps(gb, cnPgr, select="first")
foL <- findOverlaps(gb, cnLgr, select="first")

regionCN.P <- values(cnPgr)$picnic.cn[foP]
regionCN.L <- values(cnLgr)$picnic.cn[foL]


makeCompactRle <- function(reg) {
    r <- Rle(reg)
    n <- length(r@lengths)
    cbind(c(1,cumsum(r@lengths[-n])),cumsum(r@lengths-1),r@values)
}

rCH <- makeCompactRle(as.character(seqnames(gb)))
rL <- makeCompactRle(regionCN.L)
rP <- makeCompactRle(regionCN.P)

md <- rowMeans(matrix(as.numeric(rCH[,1:2]),ncol=2))

pdf("Figures/Figure_1a_SNP6_genome_wide.pdf",w=12,h=6)
par(mai=c(1.02,1.02,.1,.1))
matplot( t(rP[,1:2]), t(matrix(rep(rP[,3],2),ncol=2)), type="l", col="gray", lty=1, lwd=10,xlab="Genome Position (Chromosome)", ylab="Inferred CN state",xaxt="none", ylim=c(0,8), xaxs="i", cex.lab=2)
matplot( t(rL[,1:2]), t(matrix(rep(rL[,3],2),ncol=2)), type="l", col="black", lty=1, lwd=5, add=TRUE )
axis(1,at=md,labels=gsub("chr","",rCH[,3]), cex.axis=.8,tick=FALSE)
abline(v=as.numeric(rCH[,2]), lty=3)
dev.off()

kk <- !is.na(regionCN.L) & !is.na(regionCN.P) & regionCN.P>1
regs <- paste("L=", regionCN.L, " P=", regionCN.P, sep="")
grps <- c("L=1 P=2","L=2 P=2","L=3 P=2","L=4 P=2","L=5 P=2")

# only take well-represented groups
k <- regs %in% grps
regs[!k] <- "other"

genome.by.CN <- aggregate( width(gb)[kk], list( regs[kk] ), FUN=sum )
g <- genome.by.CN$x/sum(as.numeric(genome.by.CN$x))

pdf("Figures/Supplementary_Figure_1_cndistro.pdf",w=6,h=6)
barplot(g, names=genome.by.CN$Group, las=2, col="lightblue", ylab="Proportion")
dev.off()

gb <- gb[k]
regs <- regs[k]

save(gb,regs,file="Rdata_created/orig_regs_gb.Rdata")


## pre-processing of MBD-seq data
## --------------------------------------
#home <- "/Users/mark/projects/cn_paper"
# this reads from internal pipeline, end result is a GRangesList of reads
#setwd("/Users/mark/wehi/data/")
#f <- read.csv("NGS_Log14.csv", stringsAsFactors=FALSE)
#k <- f$Sample_Name %in% c("LNCaP","PrEC","SssI_Elution_6") & f$Experiment_Type=="MBD2IP"
#f <- f[k,]

#fn <- f$RdataGR
#names(fn) <- f$Sample_ID

#grl <- sapply(fn, function(u) { cat("."); load(u); rs })
#grl <- lapply(grl, updateObject)
#grl <- GRangesList(grl)

#setwd(home)
#save(grl, file="MBD_grl.Rdata")

# count read densities in each bin genome-wide
load("Rdata_provided/MBD_grl.Rdata")
counts <- annotationBlocksCounts(grl, gb, seq.len=300)
save(counts, file="Rdata_created/counts.Rdata")



# process the gDNA-seq into inferred CN state
# -----------------------------------------
#load("Rdata_provided/LNCAP_PREC_inputs.Rdata")
#library(BSgenome.Hsapiens36bp.UCSC.hg18mappability)
#gc.par <- GCAdjustParams(genome = Hsapiens, mappability = Hsapiens36bp,
#                         min.mappability = 50, n.bins = 20, min.bin.size = 10,
#                         poly.degree = 4, ploidy = c(2, 4))

#gbCN3 <- genomeBlocks(seqlengths(Hsapiens)[1:24], width=30000)
#counts3 <- annotationBlocksCounts(grl.inputs, gbCN3, seq.len=300)
#acn3 <- absoluteCN(gbCN3, counts3, gc.par)
#fo3 <- findOverlaps(gbCN3, gb, select="first")
#save(acn3,gbCN3,counts3,fo3,file="Rdata_created/gDNAseq_CN_LNCAP_PREC_30kb.Rdata")

load("Rdata_provided/gDNAseq_CN_LNCAP_PREC_30kb.Rdata")

cums <- cumsum(seqnames(acn3@windows)@lengths)

png("Figures/Supplementary_Figure_2_gDNA_seq_coloured_by_PICNIC.png",h=1000,w=2000)
par(mfrow=c(2,2))
plot( acn3@unadj.CN[,1], ylim=c(.5,7.5), pch=19, cex=.2, col=regionCN.P[fo3], 
      main="Unadjusted read counts (PrEC)", ylab="CN", xaxt="none", cex.main=2, xlab="")
abline(v=cums); abline(h=1:7)
plot( acn3@unadj.CN[,2], ylim=c(.5,7.5), pch=19, cex=.2, col=regionCN.L[fo3], 
      main="Unadjusted read counts (LNCaP)", ylab="CN", xaxt="none", cex.main=2, xlab="" )
abline(v=cums); abline(h=1:7)
plot( acn3@adj.CN[,1], ylim=c(.5,7.5), pch=19, cex=.2, col=regionCN.P[fo3], 
	  main="Adjusted read counts (PrEC)", ylab="CN", xaxt="none", cex.main=2, xlab="" )
abline(v=cums); abline(h=1:7)
plot( acn3@adj.CN[,2], ylim=c(.5,7.5), pch=19, cex=.2, col=regionCN.L[fo3], 
      main="Adjusted read counts (LNCaP)", ylab="CN", xaxt="none", cex.main=2, xlab="" )
abline(v=cums); abline(h=1:7)
dev.off()





## pre-processing of input-seq (gDNA-seq) data
## --------------------------------------
# this reads from internal pipeline, end result is a GRangesList of reads
# raw data is available from <DO_THIS>
#library(GenomicRanges)
#library(Repitools)
#library(BSgenome.Hsapiens.UCSC.hg18)

#home <- "/Users/mark/projects/cn_paper/infer_cn"

#setwd("/Users/mark/wehi/data/")
#f <- read.csv("NGS_Log14.csv", stringsAsFactors=FALSE)
#k <- f$Sample_Name %in% c("LNCaP","PrEC") & f$Experiment_Type=="INPUT"
#f <- f[k,]
#fn <- f$RdataGR
#names(fn) <- f$Sample_ID

#grl.inputs <- sapply(fn, function(u) { cat("."); load(u); rs })
#grl.inputs <- lapply(grl, updateObject)
#grl.inputs <- GRangesList(grl.inputs)

#setwd(home)
#save(grl, file="LNCAP_PREC_inputs.Rdata")



## pre-processing of Illumina 450k arrays
## --------------------------------------
# this reads from internal data store, end result is a GRanges of 450k array data for LNCaP/PrEC
# raw data is available from <DO_THIS>
#library(minfi)
#library(IlluminaHumanMethylation450kmanifest)

#slide.folder <- "/Users/Shared/data/450k/6264509024"
#files.table <- read.450k.sheet(slide.folder)
#files.table <-files.table[files.table[, "Sample_Group"] %in% c("PrEC", "LNCaP"), ]

#RG.raw <- read.450k.exp(base = slide.folder, targets = files.table)
#methyl.norm <- preprocessIllumina(RG.raw, bg.correct = TRUE, normalize = "controls")
#beta.table <- getBeta(methyl.norm)
#colnames(beta.table) <- files.table[, "Sample_Name"]

#meth.array.anno <- file.path(slide.folder,"..", "HumanMethylation450_15017482_v.1.2.csv")
#meth.array.anno <- read.csv(meth.array.anno, skip = 7, stringsAsFactors = FALSE)
#mapping <- match(rownames(beta.table), meth.array.anno[, "Name"])
#meth.array.anno <- meth.array.anno[mapping, ]

#k <- !is.na(as.numeric(meth.array.anno$Coordinate_36))

#beta.table <- beta.table[k,]
#meth.array.anno <- meth.array.anno[k,]

#library(GenomicRanges)
#mGR <- GRanges(seqnames=paste("chr",meth.array.anno$Chromosome_36,sep=""), 
#IRanges(start=as.numeric(meth.array.anno$Coordinate_36),width=1))

#values(mGR)$mL <- rowMeans(beta.table[,1:3])
#values(mGR)$mP <- rowMeans(beta.table[,4:6])

#save(mGR, file="illumina_450k_lncap_prec.Rdata")



# organize the H3K4me3-seq data
## --------------------------------------
# raw data is available from <DO_THIS>
#rootdir <- "/Users/Shared/data/wehi"
#fn <- dir(rootdir,"[LP].*H3K4.*Rdata",full=TRUE)
#names(fn) <- gsub(paste(rootdir,"/",sep=""),"",fn)
#names(fn) <- gsub("_GR.Rdata","",names(fn))

#grl <- sapply(fn, function(u) { cat("."); load(u); rs })
#grl <- lapply(grl, updateObject)
#grl <- GRangesList(grl)

#home <- "/Users/mark/projects/cn_paper"
#setwd(home)
#save(grl, file="H3K4me3_grl.Rdata")


# organize the H3K27me3-seq data
## --------------------------------------
# raw data is available from <DO_THIS>
#rootdir <- "/Users/Shared/data/wehi"
#fn <- dir(rootdir,"[LP].*H3K27.*Rdata",full=TRUE)
#names(fn) <- gsub(paste(rootdir,"/",sep=""),"",fn)
#names(fn) <- gsub("_GR.Rdata","",names(fn))

#grl <- sapply(fn, function(u) { cat("."); load(u); rs })
#grl <- lapply(grl, updateObject)
#grl <- GRangesList(grl)

#home <- "/Users/mark/projects/cn_paper"
#setwd(home)
#save(grl, file="H3K27me3_grl.Rdata")



