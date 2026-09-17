References:
# edgeR User's Guide: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# ggplot2 references: https://ggplot2.tidyverse.org/reference/
# ggrepel:https://ggrepel.slowkow.com/

#Load edgeR
library (edgeR)

#Load pseudobulk counts and sample metadata
pb <- readRDS("DIABETES_PSEUDOBULK/GSE241142_endothelial_pseudobulk_counts.rds")

sample_metadata <- utils::read.csv("DIABETES_PSEUDOBULK/GSE241142_endothelial_sample_metadata.csv")

#Verify dimensions and metadata
dim(pb)
sample_metadata

#Define experimental conditions, WT is used as the reference condition
condition <- factor(sample_metadata$condition, levels=c("WT", "dbdb"))
levels(condition)
table(condition)

#Create edgeR DEGList
dge_raw <- DGEList(counts=pb, group=condition)
dim(dge_raw)

#Filter low expressed genes
keep <- filterByExpr(dge_raw, group=condition)
table(keep)
dge <- dge_raw[keep,,keep.lib.sizes=FALSE]
dim(dge)

#TMM normalisation
dge <- calcNormFactors(dge)
dge$samples

#Construct design matrix
design <- model.matrix(~condition)
design
dge<- estimateDisp(dge, design)
dge$common.dispersion

#Fit quasi-likelihood model
fit <- glmQLFit(dge, design)
fit

#Test db/db versus WT
qlf <- glmQLFTest(fit, coef= "conditiondbdb")
topTags(qlf, n=5)

#Extract complete DE results
diabetes_de_all <- topTags(qlf, n=Inf, sort.by= "PValue")$table
dim(diabetes_de_all)
colnames(diabetes_de_all)

#Identify genes significant at FDR<0.05
significant_genes <- diabetes_de_all[diabetes_de_all$FDR<0.05,]
dim(significant_genes)
significant_genes

#ENSMUSG00000024222 corresponds to Fkbp5, gene symbol annotated manually from the Ensembl gene identifier
significant_genes$gene_symbol <- "Fkbp5"
significant_genes

#Construct volcano plot data preparation
volcano_data <- diabetes_de_all

#Retain Ensembl gene identifiers
volcano_data$gene <- rownames(volcano_data)

#Calculate -log10 for volcano plot y-axis
volcano_data$neg_log10_p <- -log10(volcano_data$PValue)

#Classify genes according to FDR significance
volcano_data$significance <- ifelse(volcano_data$FDR<0.05, "FDR<0.05", "Not significant")
table(volcano_data$significance)

#Labelling the significant Fkbp5 gene
volcano_data$label <- NA
volcano_data$label[volcano_data$gene == "ENSMUSG00000024222"] <- "Fkbp5"

#Check Fkbp5 values
volcano_data[volcano_data$gene == "ENSMUSG00000024222",]

Volcano plot generation
#Load the libraries 
library(ggplot2)
library(ggrepel)

#Calculate -log10 FDR for the y-axis
> volcano_data$neg_log10_FDR <- -log10(volcano_data$FDR)

#Create Volcano plot
volcano_plot <- ggplot(volcano_data,aes(x = logFC, y = neg_log10_FDR)
    
+ geom_point(aes(color = significance),
size = 1.5,
alpha = 0.6) 
+     
geom_hline(yintercept = -log10(0.05), linetype = "dashed") 
+     

 geom_vline(
xintercept = 0,
linetype = "dashed")
+     

 geom_text_repel(data = subset(volcano_data, !is.na(label)),
aes(label = label),
show.legend = FALSE,
fontface = "italic") 
+     

scale_color_manual(
values = c(
"Not significant" = "grey",
"FDR<0.05" = "red")) 
+                     
 
labs(
 x = expression(log[2]~"fold change (db/db versus WT)"),
 y = expression(-log[10]~"(FDR)")
+
  
theme_classic() +
theme(  legend.position = "none" )
