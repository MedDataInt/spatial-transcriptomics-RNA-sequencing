library(Seurat)
library(tidyverse)
library(scCustomize)  
library(ggplot2)
library(dplyr)
library(patchwork)
set.seed(2025)
setwd("/share/studies/Dermatology_Data/Data_Analysis/HS_Spatial_Transcriptomics/Jie_DataAnalysis")

object1 <- Load10X_Spatial(data.dir = "/share/studies/Dermatology_Data/Data_Analysis/HS_Spatial_Transcriptomics/Batch1_analysis/D1_1_PSX47-L5/outs", bin.size =8) # bin.size = c(2,8,16)
object2 <- Load10X_Spatial(data.dir = "/share/studies/Dermatology_Data/Data_Analysis/HS_Spatial_Transcriptomics/Batch1_analysis/D1_2_PSX41-PL2_PSX47-PL2/outs", bin.size =8) # bin.size = c(2,8,16)
object3 <- Load10X_Spatial(data.dir = "/share/studies/Dermatology_Data/Data_Analysis/HS_Spatial_Transcriptomics/Batch1_analysis/A1_1_PSX41-L9/outs", bin.size =8) # bin.size = c(2,8,16)
object4 <- Load10X_Spatial(data.dir = "/share/studies/Dermatology_Data/Data_Analysis/HS_Spatial_Transcriptomics/Batch1_analysis/A1_2_PSX62-L/outs", bin.size =8) # bin.size = c(2,8,16)

# Seurat workflow 
object1 <- NormalizeData(object1)
object1 <- FindVariableFeatures(object1)
object1 <- ScaleData(object1)


object2 <- NormalizeData(object2)
object2 <- FindVariableFeatures(object2)
object2 <- ScaleData(object2)


object3 <- NormalizeData(object3)
object3 <- FindVariableFeatures(object3)
object3 <- ScaleData(object3)

object4 <- NormalizeData(object4)
object4 <- FindVariableFeatures(object4)
object4 <- ScaleData(object4)

obj1 <- object1
obj2 <- object2
obj3 <- object3
obj4 <- object4

merged_object <- merge(object1, y = c(object2, object3, object4))


obj1 <-SketchData(obj1, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
obj2 <-SketchData(obj2, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
obj3 <-SketchData(obj3, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
obj4 <-SketchData(obj4, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")

DefaultAssay(obj1) <- "sketch"
DefaultAssay(obj2) <- "sketch"
DefaultAssay(obj3) <- "sketch"
DefaultAssay(obj4) <- "sketch"

length(rownames(obj1)) 
length(rownames(obj2)) 
length(rownames(obj3)) 
length(rownames(obj4))

#object1 <- RenameCells(object1, add.cell.id = "D1_1")
#object2 <- RenameCells(object2, add.cell.id = "D1_2")
#object3 <- RenameCells(object3, add.cell.id = "A1_1")
#object4 <- RenameCells(object4, add.cell.id = "A1_2")

obj1 <- FindVariableFeatures(obj1)
obj1 <- ScaleData(obj1)
obj1 <- RunPCA(obj1)
obj1 <- FindNeighbors(obj1, dims = 1:50)
obj1 <- FindClusters(obj1, resolution = 0.3)
obj1 <- RunUMAP(obj1, return.model = T, dims = 1:50)

obj2 <- FindVariableFeatures(obj2)
obj2 <- ScaleData(obj2)
obj2 <- RunPCA(obj2)
obj2 <- FindNeighbors(obj2, dims = 1:50)
obj2 <- FindClusters(obj2, resolution = 0.3)
obj2 <- RunUMAP(obj2, return.model = T, dims = 1:50)

obj3 <- FindVariableFeatures(obj3)
obj3 <- ScaleData(obj3)
obj3 <- RunPCA(obj3)
obj3 <- FindNeighbors(obj3, dims = 1:50)
obj3 <- FindClusters(obj3, resolution = 0.3)
obj3 <- RunUMAP(obj3, return.model = T, dims = 1:50)


obj4 <- FindVariableFeatures(obj4)
obj4 <- ScaleData(obj4)
obj4 <- RunPCA(obj4)
obj4 <- FindNeighbors(obj4, dims = 1:50)
obj4 <- FindClusters(obj4, resolution = 0.3)
obj4 <- RunUMAP(obj4, return.model = T, dims = 1:50)

saveRDS(obj1, file = "D1_1.rds")
saveRDS(obj2, file = "D1_2.rds")
saveRDS(obj3, file = "A1_1.rds")
saveRDS(obj4, file = "A1_2.rds")

#View(object@meta.data)
obj1 <- ProjectData(obj1, assay = "Spatial.008um", full.reduction = "full.pca", sketched.assay = "sketch",sketched.reduction = "pca", umap.model = "umap", dims =1:50, refdata = list(seurat_clusters.full = "seurat_clusters"))
obj2 <- ProjectData(obj2, assay = "Spatial.008um", full.reduction = "full.pca", sketched.assay = "sketch",sketched.reduction = "pca", umap.model = "umap", dims =1:50, refdata = list(seurat_clusters.full = "seurat_clusters"))
obj3 <- ProjectData(obj3, assay = "Spatial.008um", full.reduction = "full.pca", sketched.assay = "sketch",sketched.reduction = "pca", umap.model = "umap", dims =1:50, refdata = list(seurat_clusters.full = "seurat_clusters"))
obj4 <- ProjectData(obj4, assay = "Spatial.008um", full.reduction = "full.pca", sketched.assay = "sketch",sketched.reduction = "pca", umap.model = "umap", dims =1:50, refdata = list(seurat_clusters.full = "seurat_clusters"))

DefaultAssay(obj1) <- "Spatial.008um"
DefaultAssay(obj2) <- "Spatial.008um"
DefaultAssay(obj3) <- "Spatial.008um"
DefaultAssay(obj4) <- "Spatial.008um"

Idents(obj1) <- obj1@meta.data$seurat_clusters.full
Idents(obj2) <- obj2@meta.data$seurat_clusters.full
Idents(obj3) <- obj3@meta.data$seurat_clusters.full
Idents(obj4) <- obj4@meta.data$seurat_clusters.full


# List of objects
objects <- list(obj1, obj2, obj3, obj4)
names <- c("D1_1", "D1_2", "A1_1", "A1_2")

# Loop through each object and generate plots
for (i in seq_along(objects)) {
  obj <- objects[[i]]
  name <- names[i]
  
  # UMAP/tSNE clustering plot
  png(paste0(name, "_cluster_all.png"), width = 10, height = 8, units = 'in', res = 600)
  DimPlot(obj, label = TRUE)
  dev.off()
  
  # Spatial clustering plot
  png(paste0(name, "_cluster_spatial_all.png"), width = 10, height = 8, units = 'in', res = 600)
  SpatialDimPlot(obj, label = TRUE, label.size = 4)
  dev.off()
}

png("D1_1_cluster_all.png", width = 10, height = 8, units = 'in', res = 600)
DimPlot(obj1, label = T)
dev.off()

png("D1_1_cluster_spatial_all-.png", width = 10, height = 8, units = 'in', res = 600)
SpatialDimPlot(obj1, label = F, pt.size.factor = 3, label.size =4)
dev.off()


png("D1_2_cluster_all.png", width = 10, height = 8, units = 'in', res = 600)
DimPlot(obj2, label = T)
dev.off()

png("D1_2_cluster_spatial_all-.png", width = 10, height = 8, units = 'in', res = 600)
SpatialDimPlot(obj2, label = F, pt.size.factor = 3, label.size =4)
dev.off()

png("A1_1_cluster_all.png", width = 10, height = 8, units = 'in', res = 600)
DimPlot(obj3, label = T)
dev.off()

png("A1_1_cluster_spatial_all-.png", width = 10, height = 8, units = 'in', res = 600)
SpatialDimPlot(obj3, label = F, pt.size.factor = 3, label.size =20)
dev.off()


png("A1_2_cluster_all.png", width = 10, height = 8, units = 'in', res = 600)
DimPlot(obj4, label = T)
dev.off()

png("A1_2_cluster_spatial_all-.png", width = 10, height = 8, units = 'in', res = 600)
SpatialDimPlot(obj4, label = F, pt.size.factor = 3, label.size =4)
dev.off()


saveRDS(obj1, file = "D1_1a.rds")
saveRDS(obj2, file = "D1_2a.rds")
saveRDS(obj3, file = "A1_1a.rds")
saveRDS(obj4, file = "A1_2a.rds")

pdf('D1_1_PECAM1.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj1, features = "PECAM1", pt.size.factor = 8)
dev.off()
pdf('D1_2_PECAM1.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj2, features = "PECAM1", pt.size.factor = 8)
dev.off()
pdf('A1_1_PECAM1.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj3, features = "PECAM1", pt.size.factor = 8)
dev.off()
pdf('A1_2_PECAM1.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj4, features = "PECAM1", pt.size.factor = 8)
dev.off()

pdf('D1_1_DCN.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj1, features = "DCN", pt.size.factor = 8)
dev.off()
pdf('D1_2_DCN.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj2, features = "DCN", pt.size.factor = 8)
dev.off()
pdf('A1_1_DCN.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj3, features = "DCN", pt.size.factor = 8)
dev.off()
pdf('A1_2_DCN.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj4, features = "DCN", pt.size.factor = 8)
dev.off()

pdf('D1_1_KLRB1.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj1, features = "KLRB1", pt.size.factor = 8)
dev.off()
pdf('D1_2_KLRB1.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj2, features = "KLRB1", pt.size.factor = 8)
dev.off()
pdf('A1_1_KLRB1.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj3, features = "KLRB1", pt.size.factor = 8)
dev.off()
pdf('A1_2_KLRB1.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj4, features = "KLRB1", pt.size.factor = 8)
dev.off()

pdf('D1_1_ACTA2.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj1, features = "ACTA2", pt.size.factor = 8)
dev.off()
pdf('D1_2_ACTA2.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj2, features = "ACTA2", pt.size.factor = 8)
dev.off()
pdf('A1_1_ACTA2.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj3, features = "ACTA2", pt.size.factor = 8)
dev.off()
pdf('A1_2_ACTA2.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj4, features = "ACTA2", pt.size.factor = 8)
dev.off()

pdf('D1_1_CD68.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj1, features = "CD68", pt.size.factor = 8)
dev.off()
pdf('D1_2_CD68.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj2, features = "CD68", pt.size.factor = 8)
dev.off()
pdf('A1_1_CD68.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj3, features = "CD68", pt.size.factor = 8)
dev.off()
pdf('A1_2_CD68.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj4, features = "CD68", pt.size.factor = 8)
dev.off()

pdf('D1_1_IGHG1.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj1, features = "IGHG1", pt.size.factor = 8)
dev.off()
pdf('D1_2_IGHG1.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj2, features = "IGHG1", pt.size.factor = 8)
dev.off()
pdf('A1_1_IGHG1.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj3, features = "IGHG1", pt.size.factor = 8)
dev.off()
pdf('A1_2_IGHG1.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj4, features = "IGHG1", pt.size.factor = 8)
dev.off()

pdf('D1_1_CD8A.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj1, features = "CD8A", pt.size.factor = 8)
dev.off()
pdf('D1_2_CD8A.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj2, features = "CD8A", pt.size.factor = 8)
dev.off()
pdf('A1_1_CD8A.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj3, features = "CD8A", pt.size.factor = 8)
dev.off()
pdf('A1_2_CD8A.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj4, features = "CD8A", pt.size.factor = 8)
dev.off()

pdf('D1_1_CD4.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj1, features = "CD4", pt.size.factor = 8)
dev.off()
pdf('D1_2_CD4.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj2, features = "CD4", pt.size.factor = 8)
dev.off()
pdf('A1_1_CD4.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj3, features = "CD4", pt.size.factor = 8)
dev.off()
pdf('A1_2_CD4.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj4, features = "CD4", pt.size.factor = 8)
dev.off()

pdf('D1_1_CD38.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj1, features = "CD38", pt.size.factor = 8)
dev.off()
pdf('D1_2_CD38.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj2, features = "CD38", pt.size.factor = 8)
dev.off()
pdf('A1_1_CD38.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj3, features = "CD38", pt.size.factor = 8)
dev.off()
pdf('A1_2_CD38.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj4, features = "CD38", pt.size.factor = 8)
dev.off()


SpatialFeaturePlot(obj, features = c("TRAV1-2", "KLRB1", "SLC4A10", "IL18RAP")

pdf('D1_1_CD79A_.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj1, features = "CD79A", pt.size.factor = 8) + 
  scale_fill_gradient(low = "white", high = "red") +  # Adjust colors as needed
  theme_void() +  # Removes extra grid lines
  theme(panel.background = element_rect(fill = "white", color = NA))  # Sets the background to white
dev.off()

pdf('D1_2_CD79A_.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj2, features = "CD79A", pt.size.factor = 8) + 
  scale_fill_gradient(low = "white", high = "red") +  # Adjust colors as needed
  theme_void() +  # Removes extra grid lines
  theme(panel.background = element_rect(fill = "white", color = NA))  # Sets the background to white
dev.off()

pdf('A1_1_CD79A_.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj3, features = "CD79A", pt.size.factor = 8) + 
  scale_fill_gradient(low = "white", high = "red") +  # Adjust colors as needed
  theme_void() +  # Removes extra grid lines
  theme(panel.background = element_rect(fill = "white", color = NA))  # Sets the background to white
dev.off()

pdf('A1_2_CD79A_.pdf', width = 16, height = 8)
SpatialFeaturePlot(obj4, features = "CD79A", pt.size.factor = 8) + 
  scale_fill_gradient(low = "white", high = "red") +  # Adjust colors as needed
  theme_void() +  # Removes extra grid lines
  theme(panel.background = element_rect(fill = "white", color = NA))  # Sets the background to white
dev.off()

# READ
# File names and object names
file_names <- c("D1_1a.rds", "D1_2a.rds", "A1_1a.rds", "A1_2a.rds")
object_names <- c("obj1", "obj2", "obj3", "obj4")

# Read files dynamically
for (i in seq_along(file_names)) {
  assign(object_names[i], readRDS(file_names[i]))
}



genelist <- c("KRT1", "KRT10", "KRT14", "DCT", "TYRP1", "PMEL", "DCD", "MUCL1", "PIP", "PECAM1", "CDH15", "CLDN5", "DCN", "CFD", "PTGDS", "ACTA2", "TAGLN", "RGS5", "CD3D", "CD3E", "TRAC", "CD68", "LYZ", "C1QC", "TPSAB1", "TPSB2", "CTSG", "CD37", "CD69", "MS4A1", "IGHG1", "IGHG2", "IGHG3")
## heatmap
pdf('D1_1-Heatmap.pdf', width = 8, height = 10)
DoHeatmap(obj1, features = genelist) +scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
pdf('D1_1_bubble-plot.pdf', width = 8.5, height = 8)
Clustered_DotPlot(seurat_object = obj1, features = genelist)
dev.off()


pdf('D1_2-Heatmap.pdf', width = 8, height = 10)
DoHeatmap(obj2, features = genelist) +scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
pdf('D1_2_bubble-plot.pdf', width = 8.5, height = 8)
Clustered_DotPlot(seurat_object = obj2, features = genelist)
dev.off()

pdf('A1_1-Heatmap.pdf', width = 8, height = 10)
DoHeatmap(obj3, features = genelist) +scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
pdf('A1_1_bubble-plot.pdf', width = 8.5, height = 8)
Clustered_DotPlot(seurat_object = obj3, features = genelist)
dev.off()

pdf('A1_2-Heatmap.pdf', width = 8, height = 10)
DoHeatmap(obj4, features = genelist) +scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
pdf('A1_4_bubble-plot.pdf', width = 8.5, height = 8)
Clustered_DotPlot(seurat_object = obj4, features = genelist)
dev.off()
