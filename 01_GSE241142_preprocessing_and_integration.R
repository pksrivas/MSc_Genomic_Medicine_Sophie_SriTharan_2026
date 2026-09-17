library(Seurat)
setwd("~/GSE241142")

#Load the merged raw Seurat object
merged_raw <- readRDS("GSE241142_merged_raw.rds")

#load gene annotation
gene_annotation <- read.csv("GSE241142_S1_all_genes.csv.gz", check.names = FALSE)

#Mitochondrial percentage calculation
mt_genes <- gene_annotation$gene_id [grepl("^mt-|^Mt-|^MT-", gene_annotation$gene_name)]

mt_genes <- intersect(mt_genes, rownames (merged_raw))

merged_raw[["percent.mt"]] <- PercentageFeatureSet(merged_raw, features= mt_genes)

#Ribosomal percentage calculation
ribo_genes <- gene_annotation$gene_id[grepl("^Rpl|^Rps|^RPL|^RPS", gene_annotation$gene_name)]

ribo_genes <- intersect(ribo_genes, rownames(merged_raw))

merged_raw[["percent.ribo"]] <- PercentageFeatureSet(merged_raw, features = ribo_genes)

#Apply QC threshold 
after_QC <- subset(merged_raw, subset= nFeature_RNA>500 & nFeature_RNA<8000 & nCount_RNA<100000 & percent.mt<2)

#Confirm the number of nuclei retained 
dim(after_QC)

#Expected
56,980 genes x 33,713 nuclei

#SCTransform normalisation and sample integration

#Split the QC-filtered dataset by biological sample
objs <- SplitOject(after_QC, split.by= "sample_id")

#Confirm that eight biological samples are present
length(objs)

#SCTransform normalisation
#Normalise each biological sample independently using SCTransform v2 before integration.
objs <-lapply(objs, function(x) {SCTransform(x, vst.flavor="v2", verbose=FALSE)})

#Select features for integration
features <- SelectIntegrationFeatures((object.list=objs, nfeatures= 2000)
#Prepare SCT-normalised samples for integration
objs <- PrepSCTIntegration(object.list = objs, anchor.features=features)

#Identify integration anchors
#Identify transcriptionally similar nuclei across biological samples using SCT-normalised expression profiles
anchors <- FindIntegrationAnchors(object.list=objs, normalisation.method="SCT", anchor.features=features, dim=1:30, k.anchor=5, k.filter = 200, k.score=30)

#Integrate biological sample
integrated <- IntegratedData(anchorset=anchors, normalisation.method= "SCT", dims=1:30)

#Scale integrated expression data
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose= FALSE)

#Principal component analysis
#Calculate 30 principal components to reduce the dimensionality of the integrated expression matrix
integrated <- RunPCA(integrated, npcs = 30, verbose=FALSE)

#Construct nearest-neighbour graph
#Use the first 30 principal components to identify transcriptionally similar nuclei

integrated <- FindNeighbors(integrated, dims=1:30, k.param=20)

#Graph based clustering
integrated <- FindClusters(integrated, resolution=0.4)


#UMAP dimensionality reduction
#Generate 2-D representation of transcriptional relationships between nuclei
integrated <- RunUMAP(integrated, dims=1:30, n.neighbors=30, min.dist=0.3, seed.use= 42)

#Log-normalise RNA assay for gene-expression visualisation
#Return to the original RNA assay and log-normalise counts for downstream marker-expression visualisation.
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated, normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE)

#Review
dim(integrated)
Assays(integrated)
Reductions(integrated)

#Save processed integrated Seurat object
SaveRDS(integrated, file="GSE241142_integrated.rds")

#Endothelial specific reclustering
DefaultAssay(endo) <- "RNA"

endo <- NormalizeData(endo, normalization.method="LogNormalize", scale.factor=10000)

endo <- FindVariableFeatures(endo, selection.method= "vst", nfeatures=2000)

endo <- ScaleData(endo)

endo<- RunPCA (endo, npcs=50)

endo< - FindNeighbors(endo, dims=1:20, k.param=20)

endo <- FindClusters(endo, resolution=0.4)

endo <- RunUMAP(endo, dims=1:20, n.neighbors=30, min.dist=0.3, seed.use=42)


                                      
