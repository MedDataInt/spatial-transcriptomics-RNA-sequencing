# https://www.youtube.com/watch?v=zKZ03mA3aig
library(Seurat)
library(tidyverse)
library(scCustomize)  

set.seed(2025)
setwd("/share/studies/Dermatology_Data/Data_Analysis/HS_Spatial_Transcriptomics/Jie_DataAnalysis")

object1 <- Load10X_Spatial(data.dir = "/share/studies/Dermatology_Data/Data_Analysis/HS_Spatial_Transcriptomics/Batch1_analysis/D1_1_PSX47-L5/outs", bin.size =8) # bin.size = c(2,8,16)
object2 <- Load10X_Spatial(data.dir = "/share/studies/Dermatology_Data/Data_Analysis/HS_Spatial_Transcriptomics/Batch1_analysis/D1_2_PSX41-PL2_PSX47-PL2/outs", bin.size =8) # bin.size = c(2,8,16)
object3 <- Load10X_Spatial(data.dir = "/share/studies/Dermatology_Data/Data_Analysis/HS_Spatial_Transcriptomics/Batch1_analysis/A1_1_PSX41-L9/outs", bin.size =8) # bin.size = c(2,8,16)
object4 <- Load10X_Spatial(data.dir = "/share/studies/Dermatology_Data/Data_Analysis/HS_Spatial_Transcriptomics/Batch1_analysis/A1_2_PSX62-L/outs", bin.size =8) # bin.size = c(2,8,16)

Assays(object)
View(object@meta.data)

png("D1_1_QCplot.png", width = 5, height = 5, units = 'in', res = 600)
VlnPlot(object1, features = c("nCount_Spatial.008um", "nFeature_Spatial.008um"), pt.size = 0, raster = FALSE)
dev.off()

png("D1_1_QCspatial.png", width = 15, height = 5, units = 'in', res = 600)
SpatialFeaturePlot(object, features = "nCount_Spatial.008um")
dev.off()

png("D1_2_QCplot.png", width = 5, height = 5, units = 'in', res = 600)
VlnPlot(object2, features = c("nCount_Spatial.008um", "nFeature_Spatial.008um"), pt.size = 0, raster = FALSE)
dev.off()
png("A1_1_QCplot.png", width = 5, height = 5, units = 'in', res = 600)
VlnPlot(object3, features = c("nCount_Spatial.008um", "nFeature_Spatial.008um"), pt.size = 0, raster = FALSE)
dev.off()
png("A1_2_QCplot.png", width = 5, height = 5, units = 'in', res = 600)
VlnPlot(object4, features = c("nCount_Spatial.008um", "nFeature_Spatial.008um"), pt.size = 0, raster = FALSE)
dev.off()


png("D1_2_QCspatial.png", width = 15, height = 5, units = 'in', res = 600)
SpatialFeaturePlot(object2, features = "nCount_Spatial.008um")
dev.off()

png("A1_1_QCspatial.png", width = 15, height = 5, units = 'in', res = 600)
SpatialFeaturePlot(object3, features = "nCount_Spatial.008um")
dev.off()

png("A1_2_QCspatial.png", width = 15, height = 5, units = 'in', res = 600)
SpatialFeaturePlot(object4, features = "nCount_Spatial.008um")
dev.off()




# Seurat workflow 
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)

# select 50 000 cells and create a new sketch assay
object <-SketchData(object, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")

DefaultAssay(object) <- "sketch"

object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object)
object <- FindNeighbors(object, dims = 1:50)
object <- FindClusters(object, resolution = 0.5)
object <- RunUMAP(object, return.model = T, dims = 1:50)

png("D1_1_cluster.png", width = 15, height = 5, units = 'in', res = 600)
DimPlot(object, label = T)
dev.off()

png("D1_1_cluster_cus.png", width = 10, height = 8, units = 'in', res = 600)
DimPlot_scCustom(object, label = FALSE) + NoLegend()
dev.off()

View(object@meta.data)
object <- ProjectData(object, assay = "Spatial.008um", 
full.reduction = "full.pca", 
sketched.assay = "sketch",
sketched.reduction = "pca",
umap.model = "umap",
dims =1:50,
refdata = list(seurat_clusters.full = "seurat_clusters"))

# visulize the full clusters based on their spatial location 
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- object@meta.data$seurat_clusters.full

png("D1_1_cluster_all.png", width = 10, height = 8, units = 'in', res = 600)
DimPlot(object, label = T)
dev.off()

png("D1_1_cluster_spatial_all.png", width = 10, height = 8, units = 'in', res = 600)
SpatialDimPlot(object, label = T, label.size =4)
dev.off()


cells <- CellsByIdentities(object, idents = c(1,2,3))
png("D1_1_cluster_spatial_123.png", width = 10, height = 8, units = 'in', res = 600)
SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")],
cols.highlight = c("#FFF000", "grey50"),
facet.highlight = T,
combine = T)
dev.off()


png("D1_1_cluster_spatial_S100A8.png", width = 10, height = 8, units = 'in', res = 600)
SpatialFeaturePlot(object, features = "S100A8")
dev.off()

pdf('D1_1_cluster_spatial_CD38_6.pdf', width = 16, height = 8)
SpatialFeaturePlot(object, features = "CD38", pt.size.factor = 6)
dev.off()
saveRDS(object, file = "D1_1.rds")

q()
