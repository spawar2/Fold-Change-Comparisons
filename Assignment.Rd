# Author: Pawar, Date, Aug 24 2021

setwd("/Users/yalegenomecenter/Desktop")

# Reading txt file

data <- read.table(TB0004_TB0006_Treatment_samples_NumReads_genes.txt, header = TRUE, sep = "", dec = ".")

# To get gene lengths first get the gtf file from, lets say, ensembl
library(GenomicFeatures)

txdb <- makeTxDbFromGFF("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",format="gtf")

# Collect the exons per gene id

exons.list.per.gene <- exonsBy(txdb,by="gene")

# For each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})

# To calculate TPM

tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

for(i in 1:nrow(data)) {       # for-loop over rows
  data2[i, ] <- tpm3(i, exonic.gene.sizes[,1])
}

Untreated <- mean(data2$Untreated_S1, data2$Untreated_S2, data2$Untreated_S3, data2$Untreated_S4)
Stress <- mean(data2$Stress_S1, data2$Stress_S2, data2$Stress_S3, data2$Stress_S4)
both1 <- mean(data2$Drug_S1, data2$Drug_S2, data2$Drug_S3, data2$Drug_S4)
both2 <- mean(data2$Stress_Drug_S1, data2$Stress_Drug_S2, data2$Stress_Drug_S3, data2$Stress_Drug_S4)
# Cells exposed to stress compared to untreated cells
Stress/Untreated
# Cells exposed to drug compared to untreated cells
both1/Untreated
# Cells exposed to both stress and drug compared to cells exposed to stress
both2/Stress
# Cells exposed to both stress and drug compared to untreated cells
both2/Untreated
# GSEA analysis
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
df = data2
original_gene_list <- df$stats
names(original_gene_list) <- df$Gene_ID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
df2 = df[df$Gene_ID %in% dedup_ids$ENSEMBL,]
df2$Y = dedup_ids$ENTREZID
kegg_gene_list <- df2$stats
names(kegg_gene_list) <- df2$Y
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 1,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
the_plot <- cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
ridgeplot(kk2) + labs(x = "enrichment distribution")
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)

 




