# set directory
setwd("/media/Data/Microarray/POCD/")

# load packages
set.seed(1234)
library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(missMethyl)
library(doParallel)
library(lmerTest)
library(ggfortify)
library(Gviz)
library(stringr)
library(DMRcate)
library(sva)
library(ggplot2)
library(viridis)
library(bacon)
library(enrichplot)
library(wateRmelon)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# Prepare data
targets <- read.csv("POCD_samplesheet_04062023_v2_4AM.csv",header=T)
targets$Basename <- paste0(targets$Sentrix_ID,'/',targets$Sentrix_ID,'_',targets$Sentrix_Position)
targets[101:200,9:17] <- targets[1:100,9:17]
rgSet <- read.metharray.exp(targets=targets)
detP <- detectionP(rgSet)

# QC raw data
mSetRaw<- preprocessRaw(rgSet)
qcRaw <- getQC(mSetRaw)
pdf("qcRaw.pdf")
plotQC(qcRaw)
dev.off()
qcReport(rgSet)
outliers <- outlyx(rgSet,plot=TRUE)
#outliers$out
# 7th sample (Subject AC031) outlier, paired to 7th and 107th row of data. Remove both rows? 
bisConv <- bscon(rgSet)
hist(bisConv,xlab=c(0, 100))

# Removing the failed sample and their paired data
rgSet <- rgSet[,c(-7,-107)]
targets <- targets[c(-7,-107),]
detP <- detP[,c(-7,-107)]

# Predict Age and Sex to find any discrepancies
mSetRaw <- preprocessRaw(rgSet)
predAge <- agep(mSetRaw, method='all')
predSex <- estimateSex(betas(mSetRaw),do_plot=TRUE)
targets <- cbind(targets,predAge,predSex)
targets[which(targets$GENDER!=targets$predicted_sex),]
# no discrepancies in sex prediction to actual sex

# Estimate cell composition
cellCounts <- estimateCellCounts.wateRmelon(rgSet,referencePlatform="IlluminaHumanMethylationEPIC")
targets <- cbind(targets,cellCounts)

# Filter probes
# Filter CpGs
## Sex chromosomes
## SNPs
## CpH sites
## Cross-reactive
## Detection P-value > 0.01
keep <- rowSums(detP < 0.01) == ncol(mSetRaw)
mSetFiltered <- mSetRaw[keep,]
keep <- !(featureNames(mSetFiltered) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")])
mSetFiltered <- mSetFiltered[keep,]
mSetFiltered <- dropLociWithSnps(mSetFiltered)
mSetFiltered <- dropMethylationLoci(mSetFiltered)
xRtvProbes <- read.csv("epicCrossReactiveProbes.csv",header=T)
keep <- !(featureNames(mSetFiltered) %in% xRtvProbes$TargetID)
mSetFiltered <- mSetFiltered[keep,]

# Normalization
mSetFiltered.dasen <- dasen(mSetFiltered)

# Check normalization effect
qualCheck <- qual(betas(mSetFiltered), betas(mSetFiltered.dasen))
plot(qualCheck[,1],qualCheck[,2])

pdf("densityPlots.pdf")
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$CogDefWk,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$CogDefWk)),text.col=brewer.pal(8,"Set1"))
densityPlot(getBeta(mSetFiltered.dasen), sampGroups=targets$CogDefWk,main="dasen", legend=FALSE)
legend("top", legend = levels(factor(targets$CogDefWk)),text.col=brewer.pal(8,"Set1"))
dev.off()

## Generate Beta- and M-values for differential testing
mVal <- getM(mSetFiltered.dasen)
bVal <- getBeta(mSetFiltered.dasen)
save(rgSet,file="rgSet.rdata")
save(mSetFiltered.dasen,file="mSetsFiltered.dasen.rdata")
rm(mSetFiltered,keep)

# Divide up the data
targets.Ctrl <- targets[which(targets$CogDefWk==0),]
targets.POCD <- targets[which(targets$CogDefWk==1),]
bVals.Ctrl <- bVals[,which(targets$CogDefWk==0)]
bVals.POCD <- bVals[,which(targets$CogDefWk==1)]
mVals.Ctrl <- mVals[,which(targets$CogDefWk==0)]
mVals.POCD <- mVals[,which(targets$CogDefWk==1)]
targets.T1 <- targets[which(targets$Timept=="A"),]
targets.T2 <- targets[which(targets$Timept=="C"),]
bVals.T1 <- bVals[,which(targets$Timept=="A")]
bVals.T2 <- bVals[,which(targets$Timept=="C")]
mVals.T1 <- mVals[,which(targets$Timept=="A")]
mVals.T2 <- mVals[,which(targets$Timept=="C")]
