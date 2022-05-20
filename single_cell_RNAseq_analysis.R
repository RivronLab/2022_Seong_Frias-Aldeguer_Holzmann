library(dplyr)
library(Seurat)
library(patchwork)

#read the matrixes generated with FeatureCounts. The matrixes were provided by Javier Frias-Aldeguer
rivron_matrix_1 <- read.table('D:/work/rivron/matrixes/IMB-NR-001_HM5V5BGXC_extracted_aligned_assigned_sorted.tsv.gz', header=TRUE, row.names = 1)
rivron_matrix_3 <- read.table('D:/work/rivron/matrixes/IMB-NR-003_HM5V5BGXC_extracted_aligned_assigned_sorted.tsv.gz', header=TRUE, row.names = 1)
rivron_matrix_5 <- read.table('D:/work/rivron/matrixes/IMB-NR-005_HM5V5BGXC_extracted_aligned_assigned_sorted.tsv.gz', header=TRUE, row.names = 1)
rivron_matrix_6 <- read.table('D:/work/rivron/matrixes/IMB-NR-006_HM5FVBGXC_extracted_aligned_assigned_sorted.tsv.gz', header=TRUE, row.names = 1)
rivron_matrix_7 <- read.table('D:/work/rivron/matrixes/IMB-NR-007_HM5FVBGXC_extracted_aligned_assigned_sorted.tsv.gz', header=TRUE, row.names = 1)
rivron_matrix_8 <- read.table('D:/work/rivron/matrixes/IMB-NR-008_HM5FVBGXC_extracted_aligned_assigned_sorted.tsv.gz', header=TRUE, row.names = 1) #read.table

#create a Seurat object while filtering the matrix (a transcript should be found in at least 3 cells + a cell must have at least 200 transcripts)
rivron_data_1 <- CreateSeuratObject(counts = rivron_matrix_1, project = "rivron", min.cells = 3, min.features = 200)
rivron_data_3 <- CreateSeuratObject(counts = rivron_matrix_3, project = "rivron", min.cells = 3, min.features = 200)
rivron_data_5 <- CreateSeuratObject(counts = rivron_matrix_5, project = "rivron", min.cells = 3, min.features = 200)
rivron_data_6 <- CreateSeuratObject(counts = rivron_matrix_6, project = "rivron", min.cells = 3, min.features = 200)
rivron_data_7 <- CreateSeuratObject(counts = rivron_matrix_7, project = "rivron", min.cells = 3, min.features = 200)
rivron_data_8 <- CreateSeuratObject(counts = rivron_matrix_8, project = "rivron", min.cells = 3, min.features = 200)

#Create a tag for each Seurat object that will be stored in the meta.data portion of the seurat Object
#This tag will be useful to select just the cells we are interested in when we plot the cells (look at group.by = genotype in the Dimplot(PCA and UMAP))
rivron_data_1@meta.data$genotype<-"E_36_hours"
rivron_data_3@meta.data$genotype<-"F_48_hours"
rivron_data_5@meta.data$genotype<-"A_0_hours"
rivron_data_6@meta.data$genotype<-"B_4_hours"
rivron_data_7@meta.data$genotype<-"C_12_hours"
rivron_data_8@meta.data$genotype<-"D_24_hours"

#merge the datasets and create a unique matrix. 
#The add.cell.ids allow to add at the beginning of the name of each cell a name
whole<- merge(x = rivron_data_1,
              y = c(rivron_data_3, rivron_data_5, rivron_data_6, rivron_data_7, rivron_data_8 ),
                   add.cell.ids = c("36hour", "48hours", "0hours", "4hours", "12hours", "24hours"),
                   project = "rivron")

#calculate percentage mitochondrial dna in the dataset and plot it
whole[["percent.mt"]] <- PercentageFeatureSet(whole, pattern = "^mt-")
VlnPlot(whole, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(whole, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(whole, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Select only the cells that have at least 200 transcripts and at maximum 2500 + percent of mitochodrial genes < 5
#Perform the normalization and select the 10 most vaiable genes in the dataset, then create a plot
whole <- subset(whole, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
whole <- NormalizeData(whole, normalization.method = "LogNormalize", scale.factor = 10000)
whole <- FindVariableFeatures(whole, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(whole), 10)
variable_genes_1 <- VariableFeaturePlot(whole)
variable_genes_2 <- LabelPoints(plot = variable_genes_1, points = top10, repel = TRUE)
variable_genes_1 + variable_genes_2

#Run a PCA analysis and plot it
all.genes <- rownames(whole)
whole<- ScaleData(whole, features = all.genes)
whole <- RunPCA(whole, features = VariableFeatures(object = whole))
print(whole[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(whole, dims = 1:2, reduction = "pca")

DimPlot(whole, group.by = "genotype", reduction = "pca", cols = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"), pt.size = 1.5) 
DimPlot(whole, group.by = "genotype", reduction = "pca", cols = c("#F8766D", "grey", "grey", "grey", "grey", "grey"), pt.size = 1.5)
DimPlot(whole, group.by = "genotype", reduction = "pca", cols = c("grey", "#B79F00", "grey", "grey", "grey", "grey"), pt.size = 1.5)
DimPlot(whole, group.by = "genotype", reduction = "pca", cols = c("grey", "grey", "#00BA38", "grey", "grey", "grey"), pt.size = 1.5)
DimPlot(whole, group.by = "genotype", reduction = "pca", cols = c("grey", "grey", "grey", "#00BFC4", "grey", "grey"), pt.size = 1.5)
DimPlot(whole, group.by = "genotype", reduction = "pca", cols = c("grey", "grey", "grey", "grey", "#619CFF", "grey"), pt.size = 1.5)
DimPlot(whole, group.by = "genotype", reduction = "pca", cols = c("grey", "grey", "grey", "grey", "grey", "#F564E3"), pt.size = 1.5)

#Plot heat maps for the Principal components
DimHeatmap(whole, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(whole, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(whole)
whole <- FindNeighbors(whole, dims = 1:10)
whole <- FindClusters(whole, resolution = 0.5)
head(Idents(whole), 5)

#Run UMAP analysis and plot it
whole <- RunUMAP(whole, dims = 1:10)
DimPlot(whole, group.by = "genotype", reduction = "umap", cols = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"), pt.size = 1.5) 
DimPlot(whole, group.by = "genotype", reduction = "umap", cols = c("#F8766D", "grey", "grey", "grey", "grey", "grey"), pt.size = 1.5)
DimPlot(whole, group.by = "genotype", reduction = "umap", cols = c("grey", "#B79F00", "grey", "grey", "grey", "grey"), pt.size = 1.5)
DimPlot(whole, group.by = "genotype", reduction = "umap", cols = c("grey", "grey", "#00BA38", "grey", "grey", "grey"), pt.size = 1.5)
DimPlot(whole, group.by = "genotype", reduction = "umap", cols = c("grey", "grey", "grey", "#00BFC4", "grey", "grey"), pt.size = 1.5)
DimPlot(whole, group.by = "genotype", reduction = "umap", cols = c("grey", "grey", "grey", "grey", "#619CFF", "grey"), pt.size = 1.5)
DimPlot(whole, group.by = "genotype", reduction = "umap", cols = c("grey", "grey", "grey", "grey", "grey", "#F564E3"), pt.size = 1.5)


#GENES
FeaturePlot(whole, c("Cdx2", "Id2", "Fgfr1", "Pcsk6", "Myc", "Tcf7l2"), pt.size = 1)
FeaturePlot(whole, c("Igf2", "Igf2r", "Ddah1", "Dkkl1", "Mmp9", "Myc", "Car4", "Npm1", "Gas7"), pt.size = 1) #should be for 0 hours
FeaturePlot(whole, c("Plac1", "Peg10", "Ndrg1", "Ascl2", "Tbbpa", "Akt1"), pt.size = 1)
FeaturePlot(whole, c("Igf2", "Aurka", "Igf1"), pt.size = 0.5)
