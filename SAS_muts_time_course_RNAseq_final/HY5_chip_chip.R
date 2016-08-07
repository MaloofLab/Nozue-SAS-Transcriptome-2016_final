####
#
# HY5 CHIP-chip
#
#####
library(Starr)
setwd("/Volumes/Data8/GSE24973_RAW")
### sample scripts
###################################################
### code chunk number 2: Reading bpmap file
###################################################
#dataPath <- system.file("extdata", package="Starr")
#bpmapChr1 <- readBpmap(file.path(dataPath, "Scerevisiae_tlg_chr1.bpmap"))
bpmap<-readBpmap("../CD_At35b_M_v04/BPMAP/At35b_MR_v04-2_TIGRv5.bpmap") # downloaded from affymatrix web page

names(bpmap)[62]<-"Chr1"
names(bpmap)[63]<-"Chr2"
names(bpmap)[64]<-"Chr3"
names(bpmap)[65]<-"Chr4"
names(bpmap)[66]<-"Chr5"

###################################################
### code chunk number 3: Data read-in
###################################################
# cels <- c(file.path(dataPath,"Rpb3_IP_chr1.cel"), file.path(dataPath,"wt_IP_chr1.cel"), 
	# file.path(dataPath,"Rpb3_IP2_chr1.cel"))
# names <- c("rpb3_1", "wt_1","rpb3_2")
# type <- c("IP", "CONTROL", "IP")
#rpb3Chr1 <- readCelFile(bpmapChr1, cels, names, type, featureData=T, log.it=T)
cels<-c("GSM613509.CEL","GSM613510.CEL","GSM613511.CEL","GSM613512.CEL","GSM613513.CEL","GSM613514.CEL")
names<-c("WL_Hy5_rep1","WL_Hy5_rep2","WL_Hy5_rep3","WL_Hy5_rep4","WL_control_rep1","WL_control_rep2")
typ<-c("IP","IP","IP","IP","CONTROL","CONTROL")
HY5 <- readCelFile(bpmap, cels, names, type, featureData=T, log.it=T)

###################################################
### code chunk number 4: ExpressionSet
###################################################
HY5


###################################################
### code chunk number 5: assayData
###################################################
#head(exprs(rpb3Chr1))
head(exprs(HY5))


###################################################
### code chunk number 6: phenoData
###################################################
pData(rpb3Chr1)
pData(HY5)

###################################################
### code chunk number 7: featureData
###################################################
# featureData(rpb3Chr1)
featureData(HY5)

# head(featureData(rpb3Chr1)$chr)
# head(featureData(rpb3Chr1)$seq)
# head(featureData(rpb3Chr1)$pos)
head(featureData(HY5)$chr)
head(featureData(HY5)$seq)
head(featureData(HY5)$pos)


###################################################
### code chunk number 8: Reconstruction of the array image (eval = FALSE)
###################################################
## plotImage(file.path(dataPath,"Rpb3_IP_chr1.cel"))
plotImage("GSM613509.CEL")

###################################################
### code chunk number 9: Reconstruction of the array image
###################################################
# jpeg(file="image.jpeg", quality=100)
# plotImage(file.path(dataPath,"Rpb3_IP_chr1.cel"))
# dev.off()


###################################################
### code chunk number 10: boxplots and density plots (eval = FALSE)
###################################################
## par(mfcol=c(1,2))
## plotDensity(rpb3Chr1, oneDevice=T, main="")
## plotBoxes(rpb3Chr1)


###################################################
### code chunk number 11: boxplots and density plots
###################################################
pdf("boxdens.pdf", height=400, width=720)
par(mfcol=c(1,2))
plotDensity(HY5, oneDevice=T, main="")
plotBoxes(HY5)
dev.off()


###################################################
### code chunk number 12: Scatterplot matrix (eval = FALSE)
###################################################
## plotScatter(rpb3Chr1, density=T, cex=0.5)


###################################################
### code chunk number 13: Scatterplot matrix
###################################################
# png("densscatter.png", height=400, width=360)
# plotScatter(rpb3Chr1, density=T, cex=0.5)
# dev.off()


###################################################
### code chunk number 14: MA plot of raw data (eval = FALSE)
###################################################
ips <- HY5$type == "IP"
## controls <- rpb3Chr1$type == "CONTROL"
 controls <- HY5$type == "CONTROL"

## plotMA(rpb3Chr1, ip=ips, control=controls)
## plotMA(rpb3Chr1, ip=ips, control=controls)
plotMA(HY5, ip=ips, control=controls)
plotMA(HY5, ip=ips, control=controls)



###################################################
### code chunk number 15: MA plot of raw data
###################################################
# png("maRaw.png", height=400, width=720)
# ips <- rpb3Chr1$type == "IP"
# controls <- rpb3Chr1$type == "CONTROL"
# plotMA(rpb3Chr1, ip=ips, control=controls)
# dev.off()

###################################################
### code chunk number 16: Sequence-specific hybridization bias (eval = FALSE)
###################################################
## par(mfcol=c(1,2))
## plotGCbias(exprs(rpb3Chr1)[,1], featureData(rpb3Chr1)$seq, main="")
## plotPosBias(exprs(rpb3Chr1)[,1], featureData(rpb3Chr1)$seq)


###################################################
### code chunk number 17: Sequence-specific hybridization bias
###################################################
# pdf("posGC1.pdf", height=400, width=720)
# par(mfcol=c(1,2))
# plotGCbias(exprs(HY5)[,1], featureData(HY5)$seq, main="")
# plotPosBias(exprs(HY5)[,1], featureData(HY5)$seq)
# dev.off()


###################################################
### code chunk number 18: Starr.Rnw:226-227
###################################################
HY5_loess <- normalize.Probes(HY5, method="loess")


###################################################
### code chunk number 19: MA-plot of the normalized data (eval = FALSE)
###################################################
## plotMA(rpb3_loess, ip=ips, control=controls)


###################################################
### code chunk number 20: MA-plot of the normalized data
###################################################
pdf("maNorm.pdf", height=400, width=720)
plotMA(HY5_loess, ip=ips, control=controls)
dev.off()


###################################################
### code chunk number 21: Calculating ratio
###################################################
# description <- c("Rpb3vsWT")
description <- c("HY5vscontrol")
# rpb3_loess_ratio <- getRatio(rpb3_loess, ips, controls, description, fkt=median, featureData=F)
HY5_loess_ratio <- getRatio(HY5_loess, ips, controls, description, fkt=median, featureData=F)


###################################################
### code chunk number 22: Sequence-specific hybridization bias (normalized data) (eval = FALSE)
###################################################
## par(mfcol=c(1,2))
## plotGCbias(exprs(rpb3_loess_ratio)[,1], featureData(rpb3_loess)$seq, main="")
## plotPosBias(exprs(rpb3_loess_ratio)[,1], featureData(rpb3_loess)$seq, ylim=c(-0.5,0.5))


###################################################
### code chunk number 23: Sequence-specific hybridization bias (normalized data)
###################################################
# png("posGC2.png", height=400, width=720)
# par(mfcol=c(1,2))
# plotGCbias(exprs(rpb3_loess_ratio)[,1], featureData(rpb3_loess)$seq, main="")
# plotPosBias(exprs(rpb3_loess_ratio)[,1], featureData(rpb3_loess)$seq, ylim=c(-1,1))
# dev.off()


###################################################
### code chunk number 24: Starr.Rnw:298-299
###################################################

probeAnno <- bpmapToProbeAnno(bpmap)


###################################################
### code chunk number 25: Remapping probes
###################################################
#newbpmap <- remap(bpmapChr1, path=dataPath, reverse_complementary=TRUE, return_bpmap=TRUE)
newbpmap <- remap(bpmap, reverse_complementary=TRUE, return_bpmap=TRUE)
# Number of nodes: 44385334
# Searching: 
# 0 % of the probes could be mapped uniquely.
# Error in split.default(1:length(matches$text), factor(matches$text)) : 
  # group length is 0 but data length > 0


###################################################
### code chunk number 26: Summary of bpmap
###################################################
str(newbpmap)


###################################################
### code chunk number 27: Summary of bpmap (eval = FALSE)
###################################################
## writeTpmap("newbpmap.tpmap", newbpmap)
## tpmap2bpmap("newbpmap.tpmap", "newbpmap.bpmap")
## 
## pA <- bpmapToProbeAnno(newbpmap)


###################################################
### code chunk number 28: Starr.Rnw:352-354
###################################################
#transcriptAnno <- read.gffAnno(file.path(dataPath, "transcriptAnno.gff"), feature="transcript")
transcriptAnno<-read.gffAnno("TAIR10_GFF3_genes.gff",feature="mRNA")
# Error: length(take) > 0 is not TRUE

filteredIDs <- filterGenes(transcriptAnno, distance_us = 0, distance_ds = 0, minLength = 1000)


###################################################
### code chunk number 29: means
###################################################
pos <- c("start", "start", "start", "region", "region","region","region", "stop","stop","stop")
upstream <- c(500, 0, 250, 0, 0, 500, 500, 500, 0, 250)
downstream <- c(0, 500, 250, 0, 500, 0, 500, 0, 500, 250)
info <- data.frame(pos=pos, upstream=upstream, downstream=downstream, stringsAsFactors=F)


###################################################
### code chunk number 30: means
###################################################
#means_rpb3 <- getMeans(rpb3_loess_ratio, probeAnnoChr1, transcriptAnno[which(transcriptAnno$name %in% filteredIDs),], info)
means_HY5 <- getMeans(HY5_loess_ratio, probeAnno, transcriptAnno[which(transcriptAnno$name %in% filteredIDs),], info)
# Error in getProfiles(eSet, probeAnno, gffAnno = anno, upstream = regions[i,  : 
  # No chromosomes in gff annotation mapped in probeAnno object.


###################################################
### code chunk number 31: correlationPlot (eval = FALSE)
###################################################
## info$cor <- sapply(means_rpb3, mean, na.rm=T)
## level <- c(1, 1, 2, 3, 4, 5, 6, 1, 1, 2)
## info$level <- level
## correlationPlot(info, labels=c("TSS", "TTS"))


###################################################
### code chunk number 32: correlationPlot
###################################################
pdf("corPlot.pdf", height=400, width=360)
info$cor <- sapply(means_rpb3, mean, na.rm=T)
level <- c(1, 1, 2, 3, 4, 5, 6, 1, 1, 2)
info$level <- level
correlationPlot(info, labels=c("TSS", "TTS"))
dev.off()


###################################################
### code chunk number 33: profileplotExampleData (eval = FALSE)
###################################################
## sampls = 100
## probes = 63
## at = (-31:31)*14
## clus = matrix(rnorm(probes*sampls,sd=1),ncol=probes)
## clus= rbind( t(t(clus)+sin(1:probes/10))+1:nrow(clus)/sampls , t(t(clus)+sin(pi/2+1:probes/10))+1:nrow(clus)/sampls )


###################################################
### code chunk number 34: profileplotExampleData (eval = FALSE)
###################################################
## labs = paste("cluster",kmeans(clus,2)$cluster)


###################################################
### code chunk number 35: profileplotExampleData (eval = FALSE)
###################################################
## par(mfrow=c(1,2))
## profileplot(clus,label=labs,main="Clustered data",colpal=c("heat","blue"),add.quartiles=T,fromto=c(0.05,0.95))


###################################################
### code chunk number 36: profileplot
###################################################
pdf("profileplot.pdf", height=400, width=720)
sampls = 100
probes = 63
at = (-31:31)*14
clus = matrix(rnorm(probes*sampls,sd=1),ncol=probes)
clus= rbind( t(t(clus)+sin(1:probes/10))+1:nrow(clus)/sampls , t(t(clus)+sin(pi/2+1:probes/10))+1:nrow(clus)/sampls )
labs = paste("cluster",kmeans(clus,2)$cluster)
par(mfrow=c(1,2))
profileplot(clus,label=labs,main="Clustered data",colpal=c("heat","blue","red","topo"),add.quartiles=T,fromto=c(0.05,0.95))
dev.off()


###################################################
### code chunk number 37: Starr.Rnw:466-471
###################################################
tssAnno <- transcriptAnno
watson <- which(tssAnno$strand == 1)
tssAnno[watson,]$end <- tssAnno[watson,]$start
crick <- which(tssAnno$strand == -1)
tssAnno[crick,]$start <- tssAnno[crick,]$end


###################################################
### code chunk number 38: Starr.Rnw:476-477
###################################################
#profile <- getProfiles(rpb3_loess_ratio, probeAnnoChr1, tssAnno, 500, 500, feature="TSS", borderNames="TSS", method="basewise")
profile <- getProfiles(HY5_loess_ratio, probeAnno, tssAnno, 500, 500, feature="TSS", borderNames="TSS", method="basewise")


###################################################
### code chunk number 39: plotProfiles (eval = FALSE)
###################################################
## clust <- rep(1, dim(tssAnno)[1])
## names(clust) <- tssAnno$name
## plotProfiles(profile, cluster=clust)


###################################################
### code chunk number 40: plotProfiles
###################################################
pdf("sumPlot.pdf", height=400, width=720)
clust <- rep(1, dim(tssAnno)[1])
names(clust) <- tssAnno$name
plotProfiles(profile, cluster=clust, type="l", lwd=2)
dev.off()


###################################################
### code chunk number 41: Starr.Rnw:515-516
###################################################
#peaks <- cmarrt.ma(rpb3_loess_ratio, probeAnnoChr1, chr=NULL, M=NULL, frag.length=300)
peaks <- cmarrt.ma(HY5_loess_ratio, probeAnno, chr=NULL, M=NULL, frag.length=300) # does not work

# > probeAnno
# A 'probeAnno' object holding the mapping between
# reporters and genomic positions.
# Chromosomes: AffxCtrlBkGrAntiGenomic:v1;gcBin03 AffxCtrlBkGrAntiGenomic:v1;gcBin04 AffxCtrlBkGrAntiGenomic:v1;gcBin05 AffxCtrlBkGrAntiGenomic:v1;gcBin06 AffxCtrlBkGrAntiGenomic:v1;gcBin07 AffxCtrlBkGrAntiGenomic:v1;gcBin08 AffxCtrlBkGrAntiGenomic:v1;gcBin09 AffxCtrlBkGrAntiGenomic:v1;gcBin10 AffxCtrlBkGrAntiGenomic:v1;gcBin11 AffxCtrlBkGrAntiGenomic:v1;gcBin12 AffxCtrlBkGrAntiGenomic:v1;gcBin13 AffxCtrlBkGrAntiGenomic:v1;gcBin14 AffxCtrlBkGrAntiGenomic:v1;gcBin15 AffxCtrlBkGrAntiGenomic:v1;gcBin16 AffxCtrlBkGrAntiGenomic:v1;gcBin17 AffxCtrlBkGrAntiGenomic:v1;gcBin18 AffxCtrlBkGrAntiGenomic:v1;gcBin19 AffxCtrlBkGrAntiGenomic:v1;gcBin20 AffxCtrlBkGrAntiGenomic:v1;gcBin21 AffxCtrlBkGrAntiGenomic:v1;gcBin22 AffxCtrlBkGrAntiGenomic:v1;gcBin23 AffxCtrlBkGrAntiGenomic:v1;gcBin24 AffxCtrlBkGrAntiGenomic:v1;gcBin25 AffxCtrlBkGrGenomic:v1;gcBin04 AffxCtrlBkGrGenomic:v1;gcBin05 AffxCtrlBkGrGenomic:v1;gcBin06 AffxCtrlBkGrGenomic:v1;gcBin07 AffxCtrlBkGrGenomic:v1;gcBin08 AffxCtrlBkGrGenomic:v1;gcBin09 AffxCtrlBkGrGenomic:v1;gcBin10 AffxCtrlBkGrGenomic:v1;gcBin11 AffxCtrlBkGrGenomic:v1;gcBin12 AffxCtrlBkGrGenomic:v1;gcBin13 AffxCtrlBkGrGenomic:v1;gcBin14 AffxCtrlBkGrGenomic:v1;gcBin15 AffxCtrlBkGrGenomic:v1;gcBin16 AffxCtrlBkGrGenomic:v1;gcBin17 AffxCtrlBkGrGenomic:v1;gcBin18 AffxCtrlBkGrGenomic:v1;gcBin19 AffxCtrlBkGrGenomic:v1;gcBin20 AffxCtrlBkGrGenomic:v1;gcBin21 AffxCtrlBkGrGenomic:v1;gcBin22 AffxCtrlBkGrGenomic:v1;gcBin23 AffxCtrlBs:v1;r2_dap AffxCtrlBs:v1;r2_lys AffxCtrlBs:v1;r2_phe AffxCtrlBs:v1;r2_thr AffxCtrlBs:v1;X_Dap AffxCtrlBs:v1;X_Lys AffxCtrlBs:v1;X_Phe AffxCtrlBs:v1;X_Thr AffxCtrlBs:v1;X_Trpn AffxCtrlEc:v1;r2_bioB AffxCtrlEc:v1;r2_bioC AffxCtrlEc:v1;r2_bioD AffxCtrlHouseKeepingAt:v1;25S AffxCtrlHouseKeepingAt:v1;5S AffxCtrlHouseKeepingAt:v1;Actin AffxCtrlP1:v1;r2_cre AffxCtrl:v1;r2_Tag At:TIGRv5;chloroplast At:TIGRv5;mitochondria Bs:Jan_2004;NC_000964.1 Chr1 Chr2 Chr3 Chr4 Chr5
# Microarray platform:  
#   Genome: 
#   > peaks <- cmarrt.ma(HY5_loess_ratio, probeAnno, chr=NULL, M=NULL, frag.length=300)
# Error in probeAnno[paste(x, "start", sep = ".")] : 
#   No mapping 'Bs:Jan_2004;NC_000964.start' in this 'probeAnno' object.
# > HY5_loess_ratio
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 3092374 features, 1 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: HY5vscontrol
# varLabels: type
# varMetadata: labelDescription
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation:  
  



###################################################
### code chunk number 42: diagnostic plots cmarrt (eval = FALSE)
###################################################
## plotcmarrt(peaks)


###################################################
### code chunk number 43: diagnostic plots cmarrt
###################################################
pdf("cmarrt.pdf", height=800, width=720)
plotcmarrt(peaks)
dev.off()


###################################################
### code chunk number 44: Starr.Rnw:552-553
###################################################
peaklist <- cmarrt.peak(peaks, alpha = 0.05, method = "BH", minrun = 4)


###################################################
### code chunk number 45: Starr.Rnw:556-557
###################################################
str(peaklist)


###################################################
### code chunk number 46: smoothing
###################################################
# rpb3_ratio_smooth <- computeRunningMedians(rpb3_loess_ratio, probeAnno=probeAnnoChr1, allChr = "chr1", winHalfSize = 80, modColumn="type")
# sampleNames(rpb3_ratio_smooth) <- paste(sampleNames(rpb3_loess_ratio),"smoothed")
# y0 <- apply(exprs(rpb3_ratio_smooth), 2, upperBoundNull)

HY5_ratio_smooth <- computeRunningMedians(HY5_loess_ratio, probeAnno=probeAnno, winHalfSize = 80, modColumn="type")
sampleNames(HY5_ratio_smooth) <- paste(sampleNames(HY5_loess_ratio),"smoothed")
y0 <- apply(exprs(HY5_ratio_smooth), 2, upperBoundNull)

###################################################
### code chunk number 47: ChIP-enriched regions
###################################################
distCutOff <- max(transcriptAnno$end - transcriptAnno$start)
chers <- findChersOnSmoothed(HY5_ratio_smooth, probeAnno=probeAnno, thresholds=y0, distCutOff=distCutOff, minProbesInRow = 10, verbose = TRUE)

###################################################
### code chunk number 48: ChIP-enriched regions
###################################################
chers <- relateChers(chers, transcriptAnno, upstream=500)


###################################################
### code chunk number 49: ChIP-enriched regions
###################################################
chersD <- as.data.frame.cherList(chers)
chersD <- chersD[which(chersD$feature != ""),]
chersD[order(chersD$maxLevel, decreasing=TRUE)[1:5],]

# list of genes (mo)
chers.features<-unlist(strsplit(chersD$features,split=" "))
chers.features.s<-unique(gsub("([[:print:]]+)(.[[:digit:]])","\\1",chers.features)) # HY5 binding site

write.csv(chers.features.s,file="HY5_binding_site_WL.csv")

###################################################
### code chunk number 50: plotCher
###################################################
# load("temp.Rdata")
plot(chers[[86]], HY5_ratio_smooth, probeAnno=probeAnno, gff=transcriptAnno, paletteName="Spectral")

###################################################
### code chunk number 51: sessionInfo
###################################################
toLatex(sessionInfo())

###################
## by Ringo
##################

