
## Droplet processing

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##  Script overview
## This script is used to process the results from CellRanger multi. 
## You can check this 10x page(https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-5p-outputs-overview-multi) for an overview of cellranger multi output
## The script main analyses include:
## 1. Identify and analyse empty droplets 
## 2. Demultiplex the droplets and assign them to donors (check this vignette from Seurat https://satijalab.org/seurat/articles/hashing_vignette)
## 3. Prepare the data for doublet identification

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Set up the script
library(Seurat)
library(DropletUtils)
library(knitr)
library(data.table)
library(dplyr)
library(magrittr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggrastr)
library(ggrepel)
library(gridExtra)


set.seed(20240321)
options(Seurat.object.assay.version = 'v5') ## this option is important as it defines the Seurat version you are using. Note that different versions might use different functions/syntax
options(future.globals.maxSize = 1e9)

outPlotDir = createDir(path=paste(vastProjectDir$plots,scriptID,sep='/'))
outFileDir = createDir(path=paste(vastProjectDir$files,scriptID,sep='/'))

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Read cellranger output
## First read barcodes/features/matrix raw and filtered files into a Seurat object and then split the ADT,HTO and GEX into different assays.
## Splitting the ADT and HTO data into separate assays is done by reading the Antibody feature info csv file containing the antibody names

tenx_raw_out_dir = 'path/to/multi/count/raw_feature_bc_matrix/' 
features_file = 'path/to/multi/count/cellranger-multi-features.csv' ## comma-separated file containing the feature info for the ADT/HTO or any other assays

features <- fread(features_file,sep=',',header=T) ## reat the csv file into a data.table
cellranger.counts = Read10X(data.dir = tenx_raw_out_dir) ## read the count matrix (this will return a list of matrices)

## Start creating the seurat object using the Gene Expression (GEX) as main assay
seurat.object = CreateSeuratObject(counts = cellranger.counts$`Gene Expression`) 

## the lines below add the remaining count matrices to the same seurat object but as extra independent assays
seurat.object[['HTO']] = CreateAssayObject(
    counts = copy(cellranger.counts$`Antibody Capture`)[rownames(counts$`Antibody Capture`) %in% features[ id %like% 'Human_']$name, ] ## this line means: keep the rownames in the count matrix that match the 'Human_' pattern in their name
)

seurat.object[['ADT']] = CreateAssayObject(
    counts = copy(cellranger.counts$`Antibody Capture`)[ ! rownames(counts$`Antibody Capture`) %in% features[ id %like% 'Human_']$name, ] ## Similarly but opposite this line means: keep the rownames in the count matrix that DONT match the 'Human_' pattern in their name (see the ! mark)
)

seurat.object$Batch = 'Add batch info to seurat object'
seurat.object$Barcode = rownames(seurat.object@meta.data) ## Add barcode info to seurat object

## check you have the same barcodes across all assays
barcodes_gex = colnames(GetAssay(object = seurat.object, assay = "RNA"))
barcodes_adt = colnames(GetAssay(object = seurat.object, assay = "ADT"))
barcodes_hto = colnames(GetAssay(object = seurat.object, assay = "HTO"))
list.barcodes <- list(barcodes_gex,barcodes_adt,barcodes_hto)
differences = unlist(lapply(1:length(list.barcodes), function(n) setdiff(list.barcodes[[n]], unlist(list.barcodes[-n])))) ## if there are no differences this will variable will be 0

## Note that in case of differences, you need to retain only the barcodes commonly found across all assays


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Demultiplexing

## Demultiplexing is done to assign cell barcodes to samples (eg, donors). Demultiplexing step is performed separately for each batch.
## This demultiplexing strategy is based on Seurat's HTODemux workflow (https://satijalab.org/seurat/articles/hashing_vignette.html#demultiplex-cells-based-on-hto-enrichment) 
## NB: Demultiplexing is done on the filtered results from cellranger 
## So to integrate the data run a quick and rough log-normalisation of the RNA assay using `NormalizeData()`
## And then separately normalise the HTO and ADT assays using the CLR method

## prior to demultiplexing you need to perform a rough data processing which include normalisation, scaling and finding variable features
## you need to do this for each assay separately (HTO and GEX)

## GEX
seurat.object <- NormalizeData(seurat.object,assay = "RNA",normalization.method = "LogNormalize")
seurat.object <- FindVariableFeatures(seurat.object, selection.method = "mean.var.plot")
seurat.object <- ScaleData(seurat.object, features = VariableFeatures(seurat.object))

## HTO
seurat.object <- NormalizeData(seurat.object, assay = "HTO", normalization.method = "CLR")

## Actual demultiplexing
seurat.object  <- HTODemux(seurat.object, assay = "HTO", positive.quantile = 0.99,kfunc='clara')

kable(table(seurat.object$HTO_classification.global),align='c')
kable(table(seurat.object$hash.ID,seurat.object$Capture),align='c')

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Inspect demultiplexing results

kable(table(seurat.object$HTO_classification.global),align='c')
kable(table(seurat.object$hash.ID,seurat.object$Capture),align='c')


## Plot the number of cells assigned to Hastag IDs 
ggplot(seurat.object@meta.data,aes(x=Batch,fill=hash.ID))+
geom_histogram(stat='count')+
geom_text(stat='count', aes(label = after_stat(count)), position = position_stack(vjust = 0.5),size=5)+
theme_classic()+theme(legend.position='bottom')

## Demultiplexing QCs include:
## 1. Number of RNA molecules in each demultiplexed sample --> distribution should be higher for doublets and lower for negatives 
## 2. Enrichment HTO expression for each demultiplexed sample --> each sample should express only 1 hto 
Idents(seurat.object) <- "HTO_maxID"
RidgePlot(seurat.object, assay = "HTO", features = rownames(seurat.object[["HTO"]]), ncol = floor(length(rownames(seurat.object[["HTO"]]))/2))
Idents(seurat.object) <- "HTO_classification.global"
VlnPlot(seurat.object, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

## Plot singlets/doubles onto an UMAP space for visualisation sample clustering
## First remove negative cells from the object and cluster the HTO data
## then perform dimensional reduction and clustering on the data

seurat.object.subset <- subset(copy(seurat.object), idents = "Negative", invert = TRUE)
DefaultAssay(seurat.object.subset) <- "HTO"
seurat.object.subset <- ScaleData(seurat.object.subset, features = rownames(seurat.object.subset), verbose = FALSE)
seurat.object.subset <- RunPCA(seurat.object.subset, features = rownames(seurat.object.subset), approx = FALSE)
seurat.object.subset <- RunUMAP(seurat.object.subset, dims = 1:10, reduction.name = "htodemux.umap")

## Visualise results --> singlets should have their own separate clusters with doublets all in between
DimPlot(seurat.object.subset,group.by = 'HTO_classification.global',split.by='Capture')
DimPlot(seurat.object.subset,group.by = 'hash.ID',label=T)+theme(legend.position='bottom')

## Lastly, remove doublets and check how donors cluster together --> there should not be any donor-specific clusters
## Again, for this step run a quick and rough clustering analysis 
## So just use the RNA data to scale and reduce the dataset (no need for this step to look at the other modalities)
seurat.object.singlet <- subset(copy(seurat.object), idents = "Singlet")
DefaultAssay(seurat.object.singlet) <- 'RNA'
seurat.object.singlet <- FindVariableFeatures(seurat.object.singlet, selection.method = "mean.var.plot")
seurat.object.singlet <- ScaleData(seurat.object.singlet, features = VariableFeatures(seurat.object.singlet)) 
seurat.object.singlet  <- RunPCA(seurat.object.singlet, features = VariableFeatures(seurat.object.singlet))
seurat.object.singlet <- FindNeighbors(seurat.object.singlet, reduction = "pca", dims = 1:10)
seurat.object.singlet <- FindClusters(seurat.object.singlet, resolution = 0.6, verbose = FALSE)
seurat.object.singlet <- RunUMAP(seurat.object.singlet, reduction = "pca", dims = 1:10)

## visualise results 
DimPlot(seuratObj_filtered.integrated.singlet, group.by = "HTO_classification")+theme(legend.position='bottom')


##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Exporting the Seurat object containing all the barcodes (singlets/doublets/negatives) --> you will filter them later 
outRDSFile = 'path/to/the/output/rds/file'
saveRDS(seurat.object,outRDSFile)

##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()

