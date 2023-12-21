path = "~/Documents/cellbiology/brca_tcga_pan_can_atlas_2018"

data_patient  = read.delim(paste(path, "data_clinical_patient.txt", sep = "/"))

data_Rnaseq = read.delim(paste(path, "data_mrna_seq_v2_rsem.txt", sep = "/"))

data_cna = read.delim(paste(path, "data_cna.txt", sep = "/"))
rowerb = which(data_cna[,1]=="ERBB2")
# inspect cna data
hist(as.numeric(data_cna[rowerb,-c(1,2)]))

rna_cna_id = which(is.element(colnames(data_Rnaseq[,-c(1,2)]), colnames(data_cna[,-c(1,2)])))
# select only the rna cases which have cna data.
rna_cna_sub = data_Rnaseq[,2+rna_cna_id]

# Cols 1 and 2 are gene names.
assay = as.matrix(rna_cna_sub)

# Build metadata.
dim(assay)
metadata = matrix(0, dim(assay)[2],1)

colnames(data_patient)[1]

pat_ids = data_patient[,1]

cna_ids = colnames(data_cna)
length(cna_ids)

for (i in 1:dim(assay)[2]){
  pat_barcode = colnames(assay)[i]
  #pat_barcode = substr(pat_barcode, 1, 12)
  #pat_barcode = gsub("\\.", "-",pat_barcode)
  idx = which(pat_barcode == cna_ids)
  if( length(idx) == 0 )
  {
    metadata[i,1] = 0
    next
  }
  
  metadata[i,1] = 1*(as.numeric(data_cna[rowerb,idx]))
}

metadata[is.na(metadata)] =0
colnames(metadata) = "Amplify"

# inspect copy number distribution
h = hist(as.numeric(metadata[,1]),plot = FALSE)
plot(h, xaxt = 'n', xlim = c(-2,2.2),ylim= c(0,800), xlab = 'Copy-number', ylab = 'Count', main = NULL)
text(h$mids, h$counts + 20, ifelse(h$counts == 0, "", h$counts))
axis(side=1, at=seq(-2, 2, 1), labels=c('-2','-1','0','1','2'))

metadata[i,1] = 1*(as.numeric(data_cna[rowerb,idx] > 0))
# inspect ERBB2 amplification distribution
h = hist(as.numeric(metadata[,1]),plot = FALSE)
plot(h, xaxt = 'n', xlim = c(0,1.2),ylim= c(0,800), xlab = NULL, ylab = 'Count', main = NULL)
text(h$mids, h$counts + 20, ifelse(h$counts == 0, "", h$counts))
axis(side=1, at=seq(0, 1, 1), labels=c('ERBB2 Not Amplified','ERBB2 Amplified'))

library(DESeq2)
# Build DESeq Object

assay[is.na(assay)] = 0  # Impute with zeros the NA
assay[assay<0] = 0

dds <- DESeqDataSetFromMatrix(countData = round(assay),
                              colData = metadata,
                              design = ~ Amplify)

# Filter
#sum(rowSums(counts(dds) >= 10) >= 12)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize

dds <- dds[keep,]

# Normalize

dds <- DESeq(dds)

#BiocManager::install("apeglm")

resultsNames(dds)
#resLFC <- lfcShrink(dds, coef="Amplify", type="apeglm")
#resLFC

# plotDispEsts(dds)

# Get Results

res <- results(dds,alpha=.05)
plotMA(res, ylim=c(-2,2))

#plotMA(resLFC, ylim=c(-2,2))
#res

# create bins using the quantile function
# qs <- c(0, quantile(res$baseMean[res$baseMean > 0], 0:7/7))
# cut the genes into the bins
# bins <- cut(res$baseMean, qs)
# rename the levels of the bins using the middle point
# levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
# calculate the ratio of $p$ values less than .01 for each bin
# ratios <- tapply(res$pvalue, bins, function(p) mean(p < .01, na.rm=TRUE))
# plot these ratios
# barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")



# -------------------------------------------#
#rownames(resLFC) = gene_name

# Significantly Differentially Expressed
# Summary
summary(res)
gene_name = data_Rnaseq[keep,1]
rownames(res) = gene_name
signif = which(res$padj<0.05)
deg = res[signif,]


# Separate them 
dup = deg[deg[,2]>0.,]
#rowname_dup = rownames(dup)

  
dup <- dup[order(dup$log2FoldChange,decreasing = TRUE),]

dup_highmean = dup[dup$baseMean > 1000,]
dup_highmean[dup_highmean$log2FoldChange > 2,]

# rank genes by the mean of log fold change and p-value
rank_by_fc_pj <- function(dup){
  df <- as.data.frame(dup)
  df['rank_fc'] = c(1: nrow(df))
  
  df <- df[order(df$padj,decreasing = FALSE),]
  df['rank_pj'] = c(1: nrow(df))
  
  df['rank_fc_pj'] = ((df['rank_pj'] + df['rank_fc']) / 2);
  
  df <- df[order(df$rank_fc_pj,decreasing = FALSE),]
  
  df['rank_fc_pj'] <- c(1: nrow(df))
  
  df <- df[, c(1, 2, 6, 7, 9)]
  return(df)
}

df_up <- rank_by_fc_pj(dup)
df$core_enrichment

ddown = deg[deg[,2]<0.,]
ddown <- ddown[order(ddown$log2FoldChange,decreasing = FALSE),]
df_down <- rank_by_fc_pj(ddown)


#------------------------Pathway Enrichment Analysis---------------------------
# For Pathway Enrichment we need Entrez IDs
entrez_all = data_Rnaseq[keep,2][signif]
entrez_up = data_Rnaseq[keep,2][signif[deg[,2]>0.]]
entrez_down = data_Rnaseq[keep,2][signif[deg[,2]<0.]]

library(clusterProfiler)

# KEGG pathway over-representation analysis
all_paths =   enrichKEGG(gene = entrez_all, organism = 'hsa', pvalueCutoff = 0.05)

#function for gseKEGG
enrich_gse <- function(deg, entrez){
  kegg_gene_list <- deg$log2FoldChange
  names(kegg_gene_list) <- entrez
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  all_2 =   gseKEGG(geneList = kegg_gene_list, organism = 'hsa', pvalueCutoff = 0.05)
  return(all_2)
}
# transform to data_frame
display_enrich_result <- function(all_paths_2)
{
  df <- as.data.frame(all_paths_2)
  return(df)
}
all_paths_2 <- NULL
all_paths_2 <- enrich_gse(deg, entrez_all)
df_2 <- display_enrich_result(all_paths_2)
#for display
df <- display_enrich_result(all_paths)
df_display <- df[,-c(3,6,7,9)]

#parse all the genes in hsa04310
coreenrich = df['hsa04310','geneID']
enrichlist = strsplit(coreenrich, "/")[1]
enrich_info = matrix(0, length(enrichlist[[1]]) , 3)

j = 1
for(i in enrichlist[[1]] )
{
  idx = which(data_Rnaseq[,2] == i)
  
  enrich_info[j,1] = i
  enrich_info[j,2] = data_Rnaseq[idx,1]
  enrich_info[j,3] = deg[data_Rnaseq[idx,1],2]
  j = j + 1

}
# sort by absolute log fold change
enrich_info_wnt <- as.data.frame(enrich_info)
enrich_info_wnt <- enrich_info_wnt[order(abs(as.numeric(enrich_info_wnt$V3)),decreasing = TRUE),]
colnames(enrich_info_wnt) <- c('ID', 'Name', 'LogFoldchange')

# ----------------------------------------PCA------------------------------------------------#

dds$Amplify <- factor(dds$Amplify, levels = c(0,1), labels = c('ERBB2 not amplified','ERBB2 amplified'))
levels(dds$Amplify) <- c('ERBB2 not amplified','ERBB2 amplified')
plotCounts(dds, gene=which.min(res$padj),intgroup = 'Amplify')


# Transform the data to visualize
rld <- vst(dds, blind=FALSE)

# Do Principal Components Analysis and Plot
plotPCA(rld, intgroup="Amplify")

hist(assay[1000,], breaks = 50)

hist(assay(rld)[1000,], breaks = 50)
