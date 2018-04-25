biocLite()
library(SummarizedExperiment)
thca <- readRDS("seTHCA.rds")

assays(thca)
assays(thca)$counts[1:5, 1:5]
countexpr <- assays(thca)$counts
dim(countexpr)
countexpr[1:5, 1:5]
head(rownames(countexpr))
head(colnames(countexpr))
summary(colSums(countexpr))
rowData(thca)
rowRanges(thca)
colData(thca)
metadata(thca)
thca$sex
table(thca$sex)
subthca <- thca[1:5, 1:3]
dim(subthca)
assays(subthca)$counts

library(edgeR)
dge <- DGEList(counts = assays(thca)$counts,  genes = as.data.frame(rowData(thca)))
names(dge)
head(dge$samples)

ord <- order(dge$sample$lib.size)
barplot(dge$sample$lib.size[ord]/1e+06, las = 1, ylab = "Millions of reads", xlab = "Samples",
col = c("red", "blue")[thca$sex[ord]])
legend("topleft", c("female", "male"), fill = c("red", "blue"), inset = 0.01)

CPM <- t(t(dge$counts)/(dge$samples$lib.size/1e+06))
dim(CPM)
CPM[1:3, 1:7]
assays(thca)$logCPM <- cpm(dge, log = TRUE, prior.count = 0.25)
assays(thca)$logCPM[1:3, 1:7]

library(geneplotter)
par(mfrow = c(1, 2), mar = c(4, 5, 1, 1))
multidensity(as.list(as.data.frame(CPM)), xlab = "CPM", legend = NULL, main = "",
cex.axis = 1.2, cex.lab = 1.5, las = 1)
multidensity(as.list(as.data.frame(assays(thca)$logCPM)), xlab = "log2 CPM", legend = NULL,
main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)

par(mfrow = c(1, 2), mar = c(4, 5, 1, 1))
multidensity(as.list(as.data.frame(assays(thca)$logCPM)), xlab = "log2 CPM", legend = NULL,
main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)
boxplot(assays(thca)$logCPM, col = "gray", xlab = "Samples", ylab = expression(log[2] *
"CPM"), cex.axis = 1.2, cex.lab = 1.5, las = 1)

avgexp <- rowMeans(assays(thca)$logCPM)
hist(avgexp, xlab = expression(log[2] * "CPM"), main = "", las = 1, col = "gray")
abline(v = 0, col = "red", lwd = 2)

# Filtering approach 1
mask <- rowMeans(assays(thca)$logCPM) > 1
thca.filt <- thca[mask, ]
dge.filt <- dge[mask, ]
dim(thca.filt)

# Filtering approach 2
cpmcutoff <- round(10/min(dge$sample$lib.size/1e+06), digits = 1)
nsamplescutoff <- min(table(thca$sex))
mask <- rowSums(cpm(dge) > cpmcutoff) >= nsamplescutoff
thca.filt <- thca[mask, ]
dge.filt <- dge[mask, ]
dim(thca.filt)

par(mar = c(4, 5, 1, 1))
h <- hist(avgexp, xlab = expression("Expression level (" * log[2] * "CPM)"), main = "",
las = 1, col = "grey", cex.axis = 1.2, cex.lab = 1.5)
x <- cut(rowMeans(assays(thca.filt)$logCPM), breaks = h$breaks)
lines(h$mids, table(x), type = "h", lwd = 10, lend = 1, col = "darkred")
legend("topright", c("All genes", "Filtered genes"), fill = c("grey", "darkred"))

# Define groups (needed???)
dge$samples$group <- thca$sex
table(dge$samples$group)
dge.filt$samples$group <- thca$sex

plotSmear(dge, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0, col = "blue", lwd = 2)

dge.filt <- calcNormFactors(dge.filt)
head(dge.filt$samples$norm.factors)
head(dge.filt$samples$lib.size * dge.filt$samples$norm.factors)

# comparing un-normalized vs normalized
par(mfrow = c(1, 2))
plotSmear(dge, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0, col = "blue", lwd = 2)
plotSmear(dge.filt, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0, col = "blue", lwd = 2)

# comparing two groups by sample (example with 2 samples)
par(mfrow=c(1, 2), mar=c(4, 5, 3, 1))
for (i in 1:2) {
A <- rowMeans(assays(thca.filt)$logCPM) ; M <- assays(thca.filt)$logCPM[, i] - A
smoothScatter(A, M, main=colnames(thca.filt)[i], las=1, cex.axis=1.2, cex.lab=1.5, cex.main=2)
abline(h=0, col="blue", lwd=2) ; lo <- lowess(M ~ A) ; lines(lo$x, lo$y, col="red", lwd=2)
}

# Multidimensional Plots

plotMDS(dge.filt, col = c("red", "blue")[as.integer(dge.filt$samples$group)], cex = 0.7)
legend("topleft", c("female", "male"), fill = c("red", "blue"), inset = 0.05, cex = 0.7)

plotMDS(dge.filt, col = c("red", "orange", "blue")[as.integer(thca$concentration)],
cex = 0.7)
legend("topleft", levels(thca$concentration), fill = c("red", "orange", "blue"),
inset = 0.05, cex = 0.7)

# Filtering
maskbad <- colnames(thca) %in% c("NA18856", "NA18501")
dim(thca)
dim(dge.filt)
thca.filt <- thca.filt[, !maskbad]
dge.filt <- dge.filt[, !maskbad]
dim(thca.filt)

# Batch Effect

dge <- DGEList(counts = assays(thca)$counts, group = thca$sex, genes = as.data.frame(rowData(thca)))
dge <- calcNormFactors(dge)
head(colData(thca))

fc <- paste0(substr(thca$flow_cell, 1, 6), substr(thca$flow_cell, 18, 19))
# Examples of 3 surrogate variables of batch effect
table(data.frame(SEX = thca$sex, FLOWCELL = fc))
table(data.frame(SEX = thca$sex, LANE = thca$lane))
table(data.frame(SEX = thca$sex, CONCENTRATION = thca$concentration))

logCPM <- cpm(dge, log = TRUE, prior.count = 3)
d <- as.dist(1 - cor(logCPM, method = "spearman"))
batch <- as.integer(thca$concentration)
sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
names(batch) <- colnames(thca)
outcome <- as.character(thca$sex)
names(outcome) <- colnames(thca)
sampleDendrogram <- dendrapply(sampleDendrogram, function(x, batch, labels) { 
	if (is.leaf(x)) { ## for every node in the dendrogram if it is a leaf node
		attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")])) ## color by batch
		attr(x, "label") <- as.vector(labels[attr(x, "label")]) ## label by outcome
	}
	x
}, batch, outcome) ## these are the second and third arguments in the functionsampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
plot(sampleDendrogram, main = "Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch))), fill = sort(unique(batch)))

plotMDS(dge, labels = thca$sex, col = batch)
legend("bottomleft", paste("Batch", unique(batch)), fill = unique(batch), inset = 0.01)

library(sva)
mod <- model.matrix(~sex + concentration, data = colData(thca))
head(mod)
mod0 <- model.matrix(~concentration, data = colData(thca))
sv <- sva(logCPM, mod, mod0)
names(sv)
pValues <- f.pvalue(logCPM, mod, mod0)
sum(p.adjust(pValues, method = "BH") < 0.05)
hist(pValues, main = "", las = 1)
modSv <- cbind(mod, sv$sv)
mod0Sv <- cbind(mod0, sv$sv)
pValuesSv <- f.pvalue(logCPM, modSv, mod0Sv)
sum(p.adjust(pValuesSv, method = "BH") < 0.05)
hist(pValuesSv, main = "", las = 1)

mod <- model.matrix(~sex, colData(thca))
combatexp <- ComBat(logCPM, batch, mod)
dim(combatexp)
d <- as.dist(1 - cor(combatexp, method = "spearman"))
sampleClustering <- hclust(d)
sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
names(batch) <- colnames(thca)
outcome <- as.character(thca$sex)
names(outcome) <- colnames(thca)
sampleDendrogram <- dendrapply(sampleDendrogram, function(x, batch, labels) {
	if (is.leaf(x)) { ## for every node in the dendrogram if it is a leaf node
		attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")])) ## color by batch
		attr(x, "label") <- as.vector(labels[attr(x, "label")]) ## label by outcome
	}
	x
}, batch, outcome) ## these are the second and third arguments in the function
plot(sampleDendrogram, main = "Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch))), fill = sort(unique(batch)))

library(limma)
qrexp <- removeBatchEffect(logCPM, batch, design = mod)
dim(qrexp)
d <- as.dist(1 - cor(qrexp, method = "spearman"))
sampleClustering <- hclust(d)
sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
names(batch) <- colnames(thca)
outcome <- as.character(thca$sex)
names(outcome) <- colnames(thca)
sampleDendrogram <- dendrapply(sampleDendrogram, function(x, batch, labels) { 
	if (is.leaf(x)) { ## for every node in the dendrogram if it is a leaf node
		attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")])) ## color by batch
		attr(x, "label") <- as.vector(labels[attr(x, "label")]) ## label by outcome
	}
	x
}, batch, outcome) ## these are the second and third arguments in the function
plot(sampleDendrogram, main = "Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch))), fill = sort(unique(batch)))

library(corpcor)
s <- fast.svd(t(scale(t(logCPM), center = TRUE, scale = TRUE)))
pcSds <- s$d
pcSds[1] <- 0
svdexp <- s$u %*% diag(pcSds) %*% t(s$v)
colnames(svdexp) <- colnames(thca)
dim(svdexp)
d <- as.dist(1 - cor(svdexp, method = "spearman"))
sampleClustering <- hclust(d)
sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
names(batch) <- colnames(thca)
outcome <- as.character(thca$sex)
names(outcome) <- colnames(thca)
sampleDendrogram <- dendrapply(sampleDendrogram, function(x, batch, labels) {
	if (is.leaf(x)) { ## for every node in the dendrogram if it is a leaf node
		attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")])) ## color by batch
		attr(x, "label") <- as.vector(labels[attr(x, "label")]) ## label by outcome
	}
	x
}, batch, outcome) ## these are the second and third arguments in the function
plot(sampleDendrogram, main = "Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch))), fill = sort(unique(batch)))

# DE analyxis

dge <- DGEList(counts = assays(thca)$counts, group = thca$sex, genes = as.data.frame(rowData(thca)))
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE, prior.count = 0.25)
# see phenotypic data
head(colData(thca), n = 4)
table(thca$sex)

table(thca$sex, thca$concentration)
table(thca$sex, sub("_HGAC_S1000", "", thca$flow_cell))
table(thca$sex, thca$lane)

# logarithmic scale
meanUnloggedExp <- rowMeans(assays(thca)$count)
sdUnloggedExp <- apply(assays(thca)$counts, 1, sd)
meanLoggedExp <- rowMeans(logCPM)
sdLoggedExp <- apply(logCPM, 1, sd)
par(mfrow=c(1, 2))
plot(meanUnloggedExp, sdUnloggedExp, pch=".", cex=4, xlab="Unlogged mean expression", ylab="SD")
lines(lowess(meanUnloggedExp, sdUnloggedExp), lwd=2, col="red")
plot(meanLoggedExp, sdLoggedExp, pch=".", cex=4, xlab="Logged mean expression", ylab="SD")
lines(lowess(meanLoggedExp, sdLoggedExp, f=0.25), lwd=2, col="red")

# fold change
head(thca$sex, n = 8)
logCPM[1:4, 1:6]
maleExp <- logCPM[1, thca$sex == "male"]
femaleExp <- logCPM[1, thca$sex == "female"]
mean(maleExp)
mean(femaleExp)
logFC <- mean(maleExp) - mean(femaleExp)
logFC
2^logFC
maleExp <- rowMeans(logCPM[, thca$sex == "male"])
femaleExp <- rowMeans(logCPM[, thca$sex == "female"])
par(mfrow = c(1, 2))
plot(maleExp, femaleExp, xlab = "Male", ylab = "Female", pch = ".", cex = 4, las = 1)
plot((femaleExp + maleExp)/2, femaleExp - maleExp, pch = ".", cex = 4, las = 1)

# significance
log2fc <- femaleExp - maleExp
ranking <- order(abs(log2fc), decreasing = TRUE)
head(data.frame(Log2FC = round(log2fc[ranking], digits = 3), FC = round(2^log2fc[ranking],
	digits = 3), `1/FC` = round(2^(-log2fc[ranking]), digits = 3), row.names = rowData(thca)$symbol[ranking
	check.names = FALSE), n = 10)
