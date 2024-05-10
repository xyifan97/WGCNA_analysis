library(WGCNA)
library(dynamicTreeCut)
library(fastcluster)


rm(list = ls())
options(stringsAsFactors = FALSE)

datExpr = read.csv(file="L:/Loquat/WGCNA/ALL_sample/F_S_DEGs_rpkm_filter1.csv")

dim(datExpr)
names(datExpr)

#=====================================================================================
#
#  Code chunk 2 data cleaning and check
#
#=====================================================================================
datExpr0 = as.data.frame(t(datExpr))

rownames(datExpr0) = names(datExpr)[-c(1:9)]#每行对应着一个组

clean_datExpr = datExpr0[-1,]

View(clean_datExpr[1:5,1:5])

datExprMa = apply(clean_datExpr,2,as.numeric)

row.names(datExprMa) <- c("F02","F04","F06","F12","F14","F16","F32","F34","F36","S02","S04","S06","S12","S14","S16","S32","S34","S36")



gsg = goodSamplesGenes(datExprMa, verbose = 3)
gsg$allOK




#====================================================================================
#
#      
#
#====================================================================================

sampleTree = hclust(dist(datExprMa),method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)



powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExprMa, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

enableWGCNAThreads()
softPower = 20
adjacency = adjacency(datExprMa, power = softPower)



net = blockwiseModules(
        datExprMa,
        power = softPower,
        maxBlockSize = 6000,
        TOMType = "unsigned", minModuleSize = 30,
        reassignThreshold = 0, mergeCutHeight = 0.25,
        numericLabels = TRUE, pamRespectsDendro = FALSE,
        saveTOMs = F, 
        verbose = 3
)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
## assign all of the gene to their corresponding module 
## hclust for the genes.
