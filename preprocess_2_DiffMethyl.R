# POCD EPICArray Analysis

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
mSetRaw <- mapToGenome(mSetRaw)
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

# Color palette
pal <- brewer.pal(8,"Set1")

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

## data exploration of filtered, normalized data 
pdf("pcaFilteredData.pdf")
#par(mfrow=c(1,4))
plotMDS(getM(mSetFiltered.dasen), top=1000, gene.selection="common",col=pal[factor(targets$CogDefWk)], dim=c(1,2))
legend("right", legend=levels(factor(targets$CogDefWk)), text.col=pal,cex=0.7, bg="white")
plotMDS(getM(mSetFiltered.dasen), top=1000, gene.selection="common",col=pal[factor(targets$CogDefWk)], dim=c(1,3))
legend("right", legend=levels(factor(targets$CogDefWk)), text.col=pal,cex=0.7, bg="white")
plotMDS(getM(mSetFiltered.dasen), top=1000, gene.selection="common",col=pal[factor(targets$CogDefWk)], dim=c(1,4))
legend("right", legend=levels(factor(targets$CogDefWk)), text.col=pal,cex=0.7, bg="white")
plotMDS(getM(mSetFiltered.dasen), top=1000, gene.selection="common",col=pal[factor(targets$CogDefWk)], dim=c(1,5))
legend("right", legend=levels(factor(targets$CogDefWk)), text.col=pal,cex=0.7, bg="white")
plotMDS(getM(mSetFiltered.dasen), top=1000, gene.selection="common",col=pal[factor(targets$CogDefWk)], dim=c(1,6))
legend("right", legend=levels(factor(targets$CogDefWk)), text.col=pal,cex=0.7, bg="white")
plotMDS(getM(mSetFiltered.dasen), top=1000, gene.selection="common",col=pal[factor(targets$CogDefWk)], dim=c(1,7))
legend("right", legend=levels(factor(targets$CogDefWk)), text.col=pal,cex=0.7, bg="white")
plotMDS(getM(mSetFiltered.dasen), top=1000, gene.selection="common",col=pal[factor(targets$CogDefWk)], dim=c(1,8))
legend("right", legend=levels(factor(targets$CogDefWk)), text.col=pal,cex=0.7, bg="white")
plotMDS(getM(mSetFiltered.dasen), top=1000, gene.selection="common",col=pal[factor(targets$CogDefWk)], dim=c(1,9))
legend("right", legend=levels(factor(targets$CogDefWk)), text.col=pal,cex=0.7, bg="white")
plotMDS(getM(mSetFiltered.dasen), top=1000, gene.selection="common",col=pal[factor(targets$CogDefWk)], dim=c(1,10))
legend("right", legend=levels(factor(targets$CogDefWk)), text.col=pal,cex=0.7, bg="white")
dev.off()

save(mVals,file="mVals.rdata")
save(bVals,file="bVals.rdata")
annEPICSub <- annEPIC[match(rownames(mVals),annEPIC$Name),c(1,2,3,4,22,24,19)]
save(annEPIC,annEPICSub,file="annEPIC.rdata")


##################
# Differential Methylation Analysis (mixed-effects models)
##################

# load in packages and data
load("bVals.rdata")
load("targets.rdata")
load("annEPIC.rdata")
targets$Group <- paste0(targets$CogDefWk,"_",targets$Timept)
library(MASS)
library(emmeans)
library(lmerTest)
library(doParallel)
library(bacon)

# parallelize analysis
clus <- makeCluster(50)
registerDoParallel(clus)
clusterEvalQ(clus,library(lmerTest))

runLME <- function(rows,targets) {
fit <- lmer(bVals.sub[i,] ~ Group + GENDER + AGE + Neu + Mono + NK + Bcell + CD8T + CD4T + (1|Subject),data=targets,REML=F)
emod <- emmeans::emmeans(fit,"Group")
cMat <- coef(pairs(emod))[,c("c.1","c.2","c.5","c.6")]
# contrast 1 is 0_A - 0_C (difference in controls after surgery)
# contrast 2 is 0_A - 1_A (difference between controls and POCD before surgery)
# contrast 5 is 0_C - 1_C (difference between controls and POCD after surgery)
# contrast 6 is 1_A - 1_C (difference in POCD after surgery)
cResults <- as.data.frame(emmeans::emmeans(fit,"Group", contr = cMat, adjust = "none")$contrasts)
return(cResults)
}

# divide data up into segments
nSegments <- round(nrow(bVals)/1000)
nRows <- 1:nrow(bVals)
x <- split(nRows, factor(sort(rank(nRows)%%nSegments)))

# run mem on beta-values (slow…)
system.time({
res <- c()
for (ii in 1:length(x)) {
bVals.sub <- bVals[x[[ii]],]
resSub <- foreach(i=1:nrow(bVals.sub), .combine=rbind) %dopar% {
runLME(bVals.sub[i,],targets)
}
res <- rbind(res,resSub)
cat(paste0("Done with chunk ",ii,"\n"))
}
})

# split up results into each contrast
res.c1 <- res[which(res[,1]=="c.1"),]
res.c2 <- res[which(res[,1]=="c.2"),]
res.c5 <- res[which(res[,1]=="c.5"),]
res.c6 <- res[which(res[,1]=="c.6"),]

# adjust column and row names
colnames(res.c1) <- c("contrast","beta","se","df","t","pval")
colnames(res.c2) <- c("contrast","beta","se","df","t","pval")
colnames(res.c5) <- c("contrast","beta","se","df","t","pval")
colnames(res.c6) <- c("contrast","beta","se","df","t","pval")
rownames(res.c1) <- rownames(bVals)
rownames(res.c2) <- rownames(bVals)
rownames(res.c5) <- rownames(bVals)
rownames(res.c6) <- rownames(bVals)

# run bacon to adjust pvalues
res_bcn.c1 <- bacon(teststatistics=NULL,effectsizes=res.c1$beta,standarderrors=res.c1$se)
res.c1$bcn_pval <- pval(res_bcn.c1)
res.c1$bcn_se <- se(res_bcn.c1)
res.c1$bcn_beta <- es(res_bcn.c1)
res_bcn.c2 <- bacon(teststatistics=NULL,effectsizes=res.c2$beta,standarderrors=res.c2$se)
res.c2$bcn_pval <- pval(res_bcn.c2)
res.c2$bcn_se <- se(res_bcn.c2)
res.c2$bcn_beta <- es(res_bcn.c2)
res_bcn.c5 <- bacon(teststatistics=NULL,effectsizes=res.c5$beta,standarderrors=res.c5$se)
res.c5$bcn_pval <- pval(res_bcn.c5)
res.c5$bcn_se <- se(res_bcn.c5)
res.c5$bcn_beta <- es(res_bcn.c5)
res_bcn.c6 <- bacon(teststatistics=NULL,effectsizes=res.c6$beta,standarderrors=res.c6$se)
res.c6$bcn_pval <- pval(res_bcn.c6)
res.c6$bcn_se <- se(res_bcn.c6)
res.c6$bcn_beta <- es(res_bcn.c6)

# add some metadata to results
annEPICSub <- annEPIC[match(rownames(res.c1),rownames(annEPIC)),c(1,2,3,4,22,24,19)]
rownames(annEPICSub) <- rownames(res.c1)
res.c1 <- cbind(res.c1,annEPICSub)
res.c2 <- cbind(res.c2,annEPICSub)
res.c5 <- cbind(res.c5,annEPICSub)
res.c6 <- cbind(res.c6,annEPICSub)

# save results
save(res.c1,res.c2,res.c5,res.c6,file="res.rdata")
