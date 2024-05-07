library(tximport)
library(rhdf5)


base_dir <- "~/Analysis/kallisto"
sample_id <- dir(file.path(base_dir))
kal_dirs = sapply(sample_id, function(id) file.path(base_dir, id, 'abundance.tsv'))
tx.exp <- tximport(kal_dirs, type = "kallisto", txOut = TRUE)
head(tx.exp$counts)

tx2gene <- data.frame(
    TXNAME = rownames(tx.exp$counts),
    GENEID = sapply(strsplit(rownames(tx.exp$counts), '\\.'), '[', 1)
)
head(tx2gene)

gene.exp <- summarizeToGene(tx.exp, tx2gene)
head(gene.exp$counts)

gene.exp <- summarizeToGene(tx.exp, tx2gene, countsFromAbundance = "scaledTPM")
head(gene.exp$counts)


library(biomaRt)
library(tidyverse)

db <- useMart("ensembl")
Hs <- useDataset("hsapiens_gene_ensembl", mart = db)
h.all.genes = rownames(gene.exp$counts)
h.all.genes = unlist(h.all.genes)
id <- getBM(attributes = c("ensembl_transcript_id", "entrezgene_id", "hgnc_symbol"),
			filter = "ensembl_transcript_id",
			values = h.all.genes,
			mart = Hs)
exp = data.frame(ensembl_transcript_id = rownames(gene.exp$counts), gene.exp$counts)



mat = merge(id, exp, by= "ensembl_transcript_id")


write.table(mat, file = 'RNA-seq-counts.tsv', sep='	')
saveRDS(id, file = "Kallisto-id.rds")


lis = unique(mat$Gene.name)
mat2 = list()

for(i in 1:length(lis)){

	M = mat[mat$Gene.name == lis[i],]
	MM = M[,-c(1:7)]
	MMM = M[which.max(rowSums(MM)),]
	mat2 = rbind(mat2, MMM)
}


write.table(mat2, file = 'Mod-matrix-counts.tsv', sep='	')






# q-value calculation
library(qvalue)

PCAq_sel <- read_csv("Analysis/HPLMoral/PCAq_sel.csv")
hist(PCAq_sel$`p-value`)

qval <- qvalue(PCAq_sel$`p-value`, pi0 = 1)
summary(qval)
hist(qval$qvalues)

mat <- data.frame(qval$pvalues, qval$qvalues, qval$lfdr)
write.table(mat, "qval.csv", sep = ",", row.names = FALSE)






# PCA
library(ggplot2)
library(stats)
library(reshape2)

data <- mat2
pca_result <- prcomp(t(data[,c(2:7)]), scale = TRUE)
pc_scores <- as.data.frame(pca_result$x)
pc_data <- data.frame(SampleID = colnames(data[,c(2:7)]), pc_scores)


ggplot(pc_data, aes(x = PC1, y = PC2, color = SampleID)) +
  geom_point(size = 3) +
  xlab(paste("PC1 (", round(pca_result$sdev[1] / sum(pca_result$sdev) * 100, 2), "%)")) +
  ylab(paste("PC2 (", round(pca_result$sdev[2] / sum(pca_result$sdev) * 100, 2), "%)")) +
  ggtitle("PCA Analysis") +
  theme_classic()




