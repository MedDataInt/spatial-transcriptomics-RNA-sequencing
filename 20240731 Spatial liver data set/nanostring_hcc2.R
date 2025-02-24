
# download data from the cosmx-smi-human-liver-ffpe-dataset :
# go to the https://nanostring.com/resources/tiledb-array-cosmx-smi-human-liver-ffpe-dataset/, and then right click on the DOWNLOAD FILE, select copy link (showing below).
wget https://smi-public.objects.liquidweb.services/LiverDataReleaseSeurat_newUMAP.RDS
mkdir CosMx_HCC
cd CosMx_HCC
ls
unzip LiverDataReleaseTileDB.zip

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(matrixStats)
library(patchwork)
library(pheatmap)
library(Seurat)
library(RColorBrewer)
library(reshape2)
setwd("/home/jiew/CosMx_HCC/")
seurat_obj <- readRDS("LiverDataReleaseSeurat_newUMAP.RDS")
print(seurat_obj)
"""
An object of class Seurat
1197 features across 793318 samples within 3 assays
Active assay: RNA (1000 features, 0 variable features)
 2 layers present: counts, data
 2 other assays present: falsecode, negprobes
 3 dimensional reductions calculated: pca, approximateumap, approximateUMAP_bySlide
"""



pdf("test.pdf", width = 8.5, height = 8)
DimPlot(seurat_obj, reduction = "approximateumap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()



pdf("test_a.pdf", width = 8.5, height = 8)
DimPlot(seurat_obj, reduction = "approximateUMAP_bySlide", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()


# Feature plots for specific genes
pdf("test_H2AZ1.pdf", width = 8.5, height = 8)
FeaturePlot(seurat_obj, features = c("H2AZ1"), reduction = "approximateumap")
dev.off()

pdf("test_a_MKI67.pdf", width = 8, height = 4)
FeaturePlot(seurat_obj, features = c("MKI67"), reduction = "approximateUMAP_bySlide")
dev.off()

# Open PDF device for saving the plot
pdf("slide_distribution.pdf", width = 8.5, height = 8)
DimPlot(seurat_obj, reduction = "approximateumap", group.by = "slide_ID_numeric", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

pdf("test_b.pdf", width = 8.5, height = 8)
DimPlot(seurat_obj)
dev.off()


pdf("PCA.pdf",width = 12, height = 8)
PCAPlot(copy_seurat, dims = c(1, 2))
dev.off()

> # Specify the slide
  > slide <- 2
>
  > # Get the unique name for the slide
  > slideName <- unique(cellCoords$Run_Tissue_name[cellCoords$slide_ID_numeric == slide])
> slideMetadata <- metadata[metadata$Run_Tissue_name == slideName &
                              +                             metadata$fov %in% c(1:5, 22:26, 43:47),]


# Assuming 'metadata' is a data frame with the required information
slideMetadata <- metadata[metadata$Run_Tissue_name == slideName &
                            metadata$fov %in% c(1:5, 22:26, 43:47), ]

# Check the structure of slideMetadata to ensure it contains the necessary columns
str(slideMetadata)

# Assuming 'slideMetadata' contains the columns 'x_FOV_px', 'y_FOV_px', and 'cellType'
pdf("cell_type.pdf", width = 12, height = 8)
ggplot(slideMetadata, aes(x = x_FOV_px, y = y_FOV_px, color = cellType)) +
  geom_point(size = 0.3, alpha = 0.5) +
  labs(title = "Cell Types in Space", x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +
  scale_color_manual(values = scales::hue_pal()(length(unique(slideMetadata$cellType)))) +
  theme(legend.position = "right")
dev.off()

head(slot(object = seurat_obj, name = "meta.data")[2:5])

pdf("test_H2AZ1_a.pdf", width = 8.5, height = 8)
SpatialFeaturePlot(seurat_obj, features = c("H2AZ1"))
dev.off()


