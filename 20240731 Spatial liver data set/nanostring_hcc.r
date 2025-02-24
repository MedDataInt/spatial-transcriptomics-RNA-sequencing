# set up my R library 
# 
# Set R_LIBS environment variable
echo 'export R_LIBS="/home/jiew/R_libs"' >> ~/.bash_profile
source ~/.bash_profile


# download data from the cosmx-smi-human-liver-ffpe-dataset :
# go to the https://nanostring.com/resources/tiledb-array-cosmx-smi-human-liver-ffpe-dataset/, and then right click on the DOWNLOAD FILE, select copy link (showing below).
wget https://smi-public.objects.liquidweb.services/LiverDataReleaseTileDB.zip
mkdir CosMx_HCC
cd CosMx_HCC
ls
unzip LiverDataReleaseTileDB.zip

cd /home/jiew/CosMx_HCC/LiverDataRelease
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
setwd("/home/jiew/CosMx_HCC/LiverDataRelease")

setwd("/home/jiew/CosMx_HCC/Data_analysis")

# install in my personal folder 
if (!require("tiledb", quietly = TRUE))
  remotes::install_github("TileDB-Inc/TileDB-R", force = TRUE, 
  ref = "0.17.0")

if (!require("tiledbsc", quietly = TRUE))
  remotes::install_github("tiledb-inc/tiledbsc", force = TRUE, 
  ref = "8157b7d54398b1f957832f37fff0b173d355530e")

library(tiledb)
library(tiledbsc)

tiledbURI <- "/home/jiew/CosMx_HCC/LiverDataRelease"

# Reading into R, can read the arrays stored in the TileDB object using the tiledbsc R package which offers accesor functions to easily read the data into memory or access only.
# read in SOMACollection
tiledb_scdataset <- tiledbsc::SOMACollection$new(uri = tiledbURI, 
                                                 verbose = FALSE)
												 
names(tiledb_scdataset$somas)
names(tiledb_scdataset$somas$RNA$members)

# Raw Counts
# batch_mode parallelizes the readin, decreasing computation time
counts <- tiledb_scdataset$somas$RNA$X$members$counts$to_matrix(batch_mode = TRUE) 
dim(counts)

# Normalized Counts
norm <- tiledb_scdataset$somas$RNA_normalized$X$members$data$to_matrix(batch_mode = TRUE)
dim(norm)

# Cell Metadata
metadata <- tiledb_scdataset$somas$RNA$obs$to_dataframe()
dim(metadata)
metadata[1:4,1:10]


# this metadata to visualize the tissue structure by specifically plotting where each cell is in physical coordinate space.
cellCoords <- tiledb_scdataset$somas$RNA$obs$to_dataframe(
  attrs = c("x_FOV_px", "y_FOV_px", "x_slide_mm", "y_slide_mm", 
            "slide_ID_numeric", "Run_Tissue_name", "fov"))
head(cellCoords)
tail(cellCoords)

"""
x_FOV_px: X-coordinate in pixels within the field of view.
y_FOV_px: Y-coordinate in pixels within the field of view.
x_slide_mm: X-coordinate on the slide in millimeters.
y_slide_mm: Y-coordinate on the slide in millimeters.
slide_ID_numeric: Numeric identifier for the slide.
Run_Tissue_name: Name of the tissue sample.
fov: Identifier for the field of view.
"""

setwd("/home/jiew/CosMx_HCC/Data_analysis")
pdf("Tissue_structure.pdf", width = 8.5, height = 8)
ggplot(cellCoords, aes(x=x_slide_mm, y=y_slide_mm))+
  geom_point(alpha = 0.05, size = 0.01)+
  facet_wrap(~Run_Tissue_name)+
  coord_equal()+
  labs(title = "Cell coordinates in XY space")
dev.off()


# Target Transcript Coordinates
"""Transcript coordinates are x-y locations of individual transcript. 
Note that one cell contains many transcripts. These are currently stored in obsm slot. 
In future releases transcript coordinates 
will be stored in the uns slot: tiledb_scdataset$somas$RNA$uns$members$transcriptCoords """

transcriptCoords <- tiledb::tiledb_array(
  tiledb_scdataset$somas$RNA$obsm$members$transcriptCoords$uri,
  return_as="data.frame")[]
head(transcriptCoords)


#
""" One can visualize this by plotting transcripts from an individual sample 
(flow cell 1, normal liver) and a specific region within the tissue (FOV 15). 
In the graph below the centroid of each cell is shown as a black point and each 
transcript is shown as a colored point."""


slide <- 1
fov <- 15

slideName <- unique(cellCoords$Run_Tissue_name[cellCoords$slide_ID_numeric == 
                                                 slide])

fovCoords <- cellCoords[cellCoords$slide_ID_numeric == slide & 
                          cellCoords$fov == fov,]
fovTranscriptCoords <- transcriptCoords[transcriptCoords$slideID == slide & 
                                          transcriptCoords$fov == fov,]

targetCounts <- table(fovTranscriptCoords$target)

targets <- names(targetCounts[which(targetCounts >= 2500 & 
                                      targetCounts <= 3000)])
fovTranscriptCoords <- fovTranscriptCoords[fovTranscriptCoords$target %in% 
                                             targets,]

pdf("Normal_tissue_fov15.pdf", width = 8.5, height = 8)
ggplot(fovCoords, aes(x=x_FOV_px, y=y_FOV_px))+
  geom_point(alpha = 0.6, size = 0.1, color = "black")+
  geom_point(data = fovTranscriptCoords, 
             mapping = aes(x=x_FOV_px, 
                           y=y_FOV_px, 
                           color = target), 
             size = 0.3, alpha = 0.2)+
  theme_bw()+
  coord_equal()+
  guides(colour = guide_legend(override.aes = list(size=1,
                                                   alpha=1)))+
  labs(color = "RNA Target", title = paste0("RNA Transcripts in\n", 
                                            slideName, "\nFOV", fov))
dev.off()


slide <- 2 # hcc
fov <- 15

slideName <- unique(cellCoords$Run_Tissue_name[cellCoords$slide_ID_numeric == 
                                                 slide])

fovCoords <- cellCoords[cellCoords$slide_ID_numeric == slide & 
                          cellCoords$fov == fov,]
fovTranscriptCoords <- transcriptCoords[transcriptCoords$slideID == slide & 
                                          transcriptCoords$fov == fov,]

targetCounts <- table(fovTranscriptCoords$target)

targets <- names(targetCounts[which(targetCounts >= 2500 & 
                                      targetCounts <= 3000)])
fovTranscriptCoords <- fovTranscriptCoords[fovTranscriptCoords$target %in% 
                                             targets,]

pdf("HCC_tissue_fov15.pdf", width = 8.5, height = 8)
ggplot(fovCoords, aes(x=x_FOV_px, y=y_FOV_px))+
  geom_point(alpha = 0.6, size = 0.1, color = "black")+
  geom_point(data = fovTranscriptCoords, 
             mapping = aes(x=x_FOV_px, 
                           y=y_FOV_px, 
                           color = target), 
             size = 0.3, alpha = 0.2)+
  theme_bw()+
  coord_equal()+
  guides(colour = guide_legend(override.aes = list(size=1,
                                                   alpha=1)))+
  labs(color = "RNA Target", title = paste0("RNA Transcripts in\n", 
                                            slideName, "\nFOV", fov))
dev.off()

#######
# Specify the slide and field of view
slide <- 2
fov <- 9

# Get the unique name for the slide
slideName <- unique(cellCoords$Run_Tissue_name[cellCoords$slide_ID_numeric == slide])

# Filter coordinates for the specified slide and FOV
fovCoords <- cellCoords[cellCoords$slide_ID_numeric == slide & cellCoords$fov == fov, ]
fovTranscriptCoords <- transcriptCoords[transcriptCoords$slideID == slide & transcriptCoords$fov == fov, ]

# Filter for the specific gene of interest (VPS72)
gene_of_interest <- "FOXP3"
fovTranscriptCoords <- fovTranscriptCoords[fovTranscriptCoords$target == gene_of_interest, ]

# Plot the data
pdf("HCC_tissue_FOXP3.pdf", width = 8.5, height = 8)
ggplot(fovCoords, aes(x=x_FOV_px, y=y_FOV_px)) +
  geom_point(alpha = 0.6, size = 0.5, color = "black") +
  geom_point(data = fovTranscriptCoords, 
             mapping = aes(x=x_FOV_px, y=y_FOV_px, color = target), 
             size = 1, alpha = 0.8) +
  theme_bw() +
  coord_equal() +
  guides(colour = guide_legend(override.aes = list(size=1, alpha=1))) +
  labs(color = "RNA Target", title = paste0("Expression of ", gene_of_interest, " in\n", slideName, "\nFOV", fov))
dev.off()

# check if a gene expression? 
unique(transcriptCoords$target)
summary(transcriptCoords[transcriptCoords$target == "ZBTB16", ])


"""
 unique(transcriptCoords$target)
   [1] "PTK2"        "APOA1"       "TPT1"        "SERPINA1"    "TNFRSF1A"
   [6] "PSAP"        "HMGN2"       "RPL32"       "APOC1"       "BST2"
  [11] "FASN"        "SERPINA3"    "GLUL"        "FGG"         "CTNNB1"
  [16] "CRYAB"       "TIMP1"       "NEAT1"       "TTR"         "APOE"
  [21] "LY75"        "GC"          "QRFPR"       "VTN"         "CSTB"
  [26] "CRP"         "MALAT1"      "XBP1"        "CLU"         "LGALS3BP"
  [31] "DUSP1"       "SLCO2B1"     "COL5A2"      "ATP5F1B"     "CD81"
  [36] "NACA"        "ADGRF5"      "RPL34"       "IFIH1"       "FAU"
  [41] "TNFSF10"     "RPL37"       "SELENOP"     "ITK"         "PPIA"
  [46] "ITGAV"       "SQSTM1"      "RPL21"       "CD63"        "BAG3"
  [51] "CFD"         "EIF5A"       "PTEN"        "FN1"         "RXRA"
  [56] "H4C3"        "KRT86"       "CCL2"        "SAT1"        "C11orf96"
  [61] "AZGP1"       "UBA52"       "INSIG1"      "ARHGDIB"     "CXCL14"
  [66] "IL23A"       "CDH1"        "CXCL12"      "MYL9"        "ETV5"
  [71] "IFI6"        "CCR1"        "HLA-A"       "FCER1G"      "PLD3"
  [76] "CD5L"        "TAGLN"       "CXCR3"       "LYVE1"       "CD14"
  [81] "S100P"       "IL17D"       "B2M"         "KRT8"        "CD68"
  [86] "SREBF1"      "PPARA"       "RPL22"       "RARRES2"     "COL9A2"
  [91] "HSPA1A"      "PTTG1"       "RSPO3"       "CTSD"        "JUNB"
  [96] "TUBB"        "NUPR1"       "DUSP6"       "PTGES3"      "MIF"
 [101] "HBA1"        "COL4A5"      "ST6GAL1"     "CXCL9"       "EFNB2"
 [106] "G0S2"        "VEGFA"       "FGF2"        "DPP4"        "MT2A"
 [111] "ITGAE"       "KRT17"       "MT1X"        "IL32"        "APP"
 [116] "TNF"         "IL16"        "GLUD1"       "PGK1"        "SAA1"
 [121] "SORBS1"      "HMGCS1"      "BGN"         "BMP1"        "TSC22D1"
 [126] "HDAC11"      "CXCR2"       "LINC01857"   "PLAC8"       "TUBB4B"
 [131] "IFNAR2"      "ICOSLG"      "IL1R1"       "TNFAIP6"     "PXDN"
 [136] "PFN1"        "CCL22"       "EPHB6"       "ATP5F1E"     "JUN"
 [141] "ANXA2"       "FAS"         "RPS4Y1"      "IL6ST"       "FGFR3"
 [146] "CALD1"       "HLA-DPB1"    "TNFRSF13B"   "CCND1"       "GSN"
 [151] "COL27A1"     "YBX3"        "IL1R2"       "IGFBP7"      "ATG10"
 [156] "DCN"         "COL3A1"      "COL6A2"      "RBPJ"        "RYR2"
 [161] "FAM30A"      "BTG1"        "AIF1"        "CUZD1"       "SOD2"
 [166] "CST7"        "IL1A"        "TNXB"        "CX3CR1"      "ITM2B"
 [171] "IFITM3"      "LAMP2"       "IL6"         "CXCL17"      "IGF2"
 [176] "MTOR"        "MFAP5"       "SH3BGRL3"    "C1QB"        "RAC1"
 [181] "BECN1"       "COL18A1"     "MYL12A"      "LDLR"        "ERBB3"
 [186] "CLEC4D"      "ACKR3"       "CYP2U1"      "INSR"        "CLOCK"
 [191] "C1QC"        "SELPLG"      "KRT18"       "CLEC2D"      "EZR"
 [196] "HSP90AA1"    "MX1"         "BMPR2"       "MMP12"       "IFIT3"
 [201] "TYMS"        "AR"          "RARRES1"     "S100A6"      "CXCL8"
 [206] "IAPP"        "CALM2"       "SRGN"        "CALM1"       "IGF2R"
 [211] "FLT1"        "CLEC4E"      "OASL"        "GPX1"        "CD3D"
 [216] "INHA"        "ADIRF"       "TPM2"        "CXCR5"       "RYK"
 [221] "VIM"         "HLA-DRB1"    "PTPRC"       "PSD3"        "CD3G"
 [226] "FGF12"       "NOTCH1"      "FOS"         "ADGRA3"      "BAX"
 [231] "OAS2"        "IGFBP3"      "AXL"         "JAK1"        "PROX1"
 [236] "NLRP2"       "RARA"        "IGFBP5"      "TNFSF13B"    "SOD1"
 [241] "THBS2"       "CCR2"        "WNT10B"      "CHI3L1"      "WNT5A"
 [246] "DDC"         "AHI1"        "TNFRSF12A"   "ITGB5"       "GPX3"
 [251] "DNMT1"       "AATK"        "CD59"        "KLRF1"       "CYTOR"
 [256] "LAIR1"       "ARF1"        "HSP90B1"     "KRAS"        "ENO1"
 [261] "CD1C"        "CSPG4"       "TLR5"        "STAT1"       "MRC2"
 [266] "SQLE"        "CASP8"       "ACP5"        "LDHA"        "HSP90AB1"
 [271] "FGR"         "CAV1"        "PDS5A"       "COL17A1"     "HLA-DQB1"
 [276] "STAT5B"      "LTBR"        "RAMP1"       "CXCL5"       "ZFP36"
 [281] "RBM47"       "CD209"       "CLEC1A"      "HEXB"        "TGFBI"
 [286] "BST1"        "MARCO"       "HGF"         "TP53"        "IL15"
 [291] "CYSTM1"      "ABL2"        "LEFTY1"      "CELSR1"      "DLL1"
 [296] "NR3C1"       "CARMN"       "IL11"        "ADGRG6"      "NR1H3"
 [301] "VHL"         "ACACB"       "ITGA8"       "CD74"        "COL15A1"
 [306] "PDCD1LG2"    "C1QA"        "EPHA4"       "CCL19"       "IGF1"
 [311] "IL2RG"       "MAF"         "LGALS9"      "SLC40A1"     "NANOG"
 [316] "TXK"         "CFLAR"       "HEY1"        "IL13RA1"     "GATA3"
 [321] "DDIT3"       "ATG12"       "PDGFD"       "ITGA1"       "ACTG2"
 [326] "ZBTB16"      "BTF3"        "NPPC"        "GPBAR1"      "VWA1"
 [331] "FPR1"        "RACK1"       "ANXA4"       "ITGB1"       "B3GNT7"
 [336] "RELA"        "TNFSF14"     "ERBB2"       "LYN"         "CHEK1"
 [341] "ADGRL4"      "BTK"         "PGR"         "THBS1"       "FCGR3A"
 [346] "NR1H2"       "SOX4"        "CASP3"       "IL4R"        "MAP1LC3B"
 [351] "ACTA2"       "ARG1"        "COL5A3"      "SCGB3A1"     "IL17RB"
 [356] "EGF"         "THSD4"       "TPM1"        "CXCL1"       "HPGDS"
 [361] "TLR2"        "SPARCL1"     "S100A9"      "MAML2"       "CASR"
 [366] "CD164"       "HBB"         "ADGRL2"      "INHBA"       "NFKBIA"
 [371] "TTN"         "CD24"        "TPI1"        "GNLY"        "MARCKSL1"
 [376] "SERPINB5"    "AGR2"        "FYN"         "PECAM1"      "FABP5"
 [381] "EGFR"        "EFNA1"       "IL1RL1"      "DST"         "CD163"
 [386] "ITGA5"       "CDKN1A"      "CLEC2B"      "VCAM1"       "NRG1"
 [391] "NLRC4"       "ACVR1"       "BMP4"        "ACKR1"       "CCL20"
 [396] "IL10RB"      "ATR"         "IL6R"        "SOX9"        "RGS2"
 [401] "COL4A2"      "LTF"         "WNT3"        "FGF9"        "SIGIRR"
 [406] "LGALS1"      "CCL3"        "DHRS2"       "PHLDA2"      "IL10RA"
 [411] "SRC"         "PTGDR2"      "ITGB4"       "SMO"         "SOSTDC1"
 [416] "NGFR"        "IL3RA"       "BRCA1"       "TNFRSF14"    "HCAR2"
 [421] "SPRY4"       "MZT2A"       "PDGFC"       "SEC23A"      "NPR3"
 [426] "CELSR2"      "HCK"         "YES1"        "STAT3"       "ITGAX"
 [431] "LY6D"        "TYROBP"      "SERPINH1"    "SYK"         "PPARG"
 [436] "EPOR"        "IGKC"        "TNFRSF4"     "ADM2"        "KRT1"
 [441] "KRT80"       "FES"         "PNOC"        "CCR7"        "FFAR3"
 [446] "CXCL16"      "IL12A"       "RELT"        "CXCL13"      "SRSF2"
 [451] "EZH2"        "SCG5"        "TOP2A"       "DUSP5"       "CCL15"
 [456] "TNFRSF1B"    "KRT20"       "CIDEA"       "CD55"        "DUSP4"
 [461] "RGCC"        "KRT16"       "NR2F2"       "PIGR"        "FGF13"
 [466] "TCF7"        "BASP1"       "FOXP3"       "TNFSF9"      "WIF1"
 [471] "COL14A1"     "TNNC1"       "TEK"         "CCL21"       "DNMT3A"
 [476] "HDAC4"       "HLA-DRA"     "IFNGR2"      "RUNX3"       "MXRA8"
 [481] "PDGFB"       "ANXA1"       "IGHD"        "ITGB2"       "IGHG1"
 [486] "WNT5B"       "CD34"        "FZD5"        "IL7"         "ESR1"
 [491] "IGHA1"       "ETV4"        "SMARCB1"     "AKT1"        "IFITM1"
 [496] "CIITA"       "CCL5"        "LGR5"        "ACVR2A"      "HSPG2"
 [501] "RXRB"        "OSMR"        "COL6A1"      "TGFB3"       "LIF"
 [506] "LAMP3"       "FKBP11"      "CSF2RA"      "NFKB1"       "TLR4"
 [511] "MAP2K1"      "SPOCK2"      "NDUFA4L2"    "MET"         "CD8B"
 [516] "TNFRSF21"    "SMAD2"       "GADD45B"     "BIRC3"       "CD58"
 [521] "CCDC80"      "CD86"        "C5AR2"       "CD19"        "KLK3"
 [526] "RB1"         "DPT"         "HSPB1"       "BID"         "KLRK1"
 [531] "CLDN4"       "CCL4"        "SMAD4"       "LINC01781"   "TNFRSF10A"
 [536] "IL20"        "ACVRL1"      "CNTFR"       "PLCG1"       "SFN"
 [541] "HSD17B2"     "FGFR2"       "IL1B"        "CD84"        "MYC"
 [546] "CCL13"       "CD38"        "PDGFRA"      "ICAM1"       "HLA-DPA1"
 [551] "PCNA"        "POU5F1"      "IL20RA"      "CD93"        "ICOS"
 [556] "BRAF"        "ABL1"        "EPHB4"       "IL22RA1"     "ELANE"
 [561] "CD274"       "FKBP5"       "CD47"        "CD40"        "BBLN"
 [566] "EFNA4"       "CD8A"        "CCR5"        "CSF3"        "EPHA3"
 [571] "SNAI2"       "AHR"         "TOX"         "CD44"        "FFAR4"
 [576] "CALM3"       "HLA-DQA1"    "ALCAM"       "CACNA1C"     "CSF1R"
 [581] "ITGA6"       "SLPI"        "TNFRSF10B"   "LPAR5"       "ATM"
 [586] "BMP5"        "FZD6"        "LINC02446"   "ADGRF1"      "NKG7"
 [591] "GSK3B"       "IFIT1"       "ESAM"        "HDAC1"       "MS4A6A"
 [596] "GPNMB"       "SNAI1"       "AQP3"        "IRF4"        "INS"
 [601] "TNFRSF9"     "OAS1"        "IL15RA"      "ACVR1B"      "NCAM1"
 [606] "BCL2"        "IL12RB1"     "CD37"        "S100A4"      "HDAC3"
 [611] "BCL2L1"      "CCR10"       "KRT15"       "LIFR"        "CSF3R"
 [616] "IFNAR1"      "FZD7"        "PF4"         "CMKLR1"      "MAPK13"
 [621] "ICAM3"       "SPP1"        "XCL1"        "RARB"        "CCL28"
 [626] "KRT19"       "OSM"         "MAPK14"      "LGALS3"      "TNFSF4"
 [631] "GAS6"        "PRSS2"       "COL1A1"      "CD80"        "IL12RB2"
 [636] "S100A2"      "MMP1"        "NDRG1"       "LTB"         "ICAM2"
 [641] "CSHL1"       "GSTP1"       "STAT5A"      "SOX2"        "STAT6"
 [646] "OLFM4"       "IL24"        "CLEC10A"     "COL6A3"      "HDAC5"
 [651] "CDH19"       "ATG5"        "FZD8"        "HIF1A"       "EFNB1"
 [656] "TIGIT"       "KRT6A"       "KLF2"        "DDX58"       "IL18R1"
 [661] "IL17B"       "MST1R"       "RARG"        "CD9"         "FOXF1"
 [666] "ACKR4"       "IL17RA"      "IL17A"       "CHEK2"       "BMP2"
 [671] "GZMK"        "KRT14"       "PRF1"        "CD79A"       "COL9A3"
 [676] "ATF3"        "ITGAM"       "CCL17"       "MSR1"        "NPR2"
 [681] "FZD1"        "PTPRCAP"     "TNFRSF18"    "NCR1"        "CD33"
 [686] "CSF2RB"      "PPARD"       "RORA"        "CD28"        "RNF43"
 [691] "LYZ"         "CLEC4A"      "TNFSF12"     "IGF1R"       "CD40LG"
 [696] "MEG3"        "GDF15"       "IFNL3"       "LAG3"        "FCGBP"
 [701] "OLR1"        "RAC2"        "ST6GALNAC3"  "CALB1"       "CD4"
 [706] "IFNGR1"      "IGFBP6"      "MMP7"        "IER3"        "ENG"
 [711] "PLAC9"       "ICA1"        "PSCA"        "KRT4"        "HCST"
 [716] "MYL4"        "EPHB2"       "SPRY2"       "CCL11"       "SMAD3"
 [721] "ANGPT2"      "H2AZ1"       "ETS1"        "ADGRE2"      "MMP2"
 [726] "LDB2"        "DDR2"        "COL21A1"     "BIRC5"       "VSIR"
 [731] "ADGRA2"      "IL7R"        "PTGES"       "UBE2C"       "CD276"
 [736] "CXCR6"       "MMP19"       "LCN2"        "S100A10"     "TCAP"
 [741] "ADGRG5"      "ISG15"       "CXCL10"      "COL5A1"      "IL10"
 [746] "NOTCH3"      "G6PD"        "DMBT1"       "PTGDS"       "IFI27"
 [751] "TNFRSF11B"   "ADGRG3"      "PTK6"        "KIT"         "FASLG"
 [756] "PDGFRB"      "PTGES2"      "MERTK"       "TM4SF1"      "IRF3"
 [761] "EFNA5"       "SEC61G"      "ANGPTL1"     "IL18"        "WNT7B"
 [766] "LMNA"        "ITGAL"       "COX4I2"      "CSF1"        "TSHZ2"
 [771] "EPHA2"       "COTL1"       "SPINK1"      "ADGRL1"      "IL33"
 [776] "FFAR2"       "TYK2"        "MPO"         "VEGFD"       "ADGRE5"
 [781] "KITLG"       "KDR"         "CTLA4"       "TGFBR2"      "EPHB3"
 [786] "PDCD1"       "CYP1B1"      "CEACAM1"     "CD70"        "RGS1"
 [791] "FGF1"        "IFI44L"      "FZD3"        "NLRP1"       "ROR1"
 [796] "COL11A1"     "MS4A4A"      "APOD"        "RAMP3"       "IL1RN"
 [801] "KRT7"        "PDGFA"       "IKZF3"       "S100A8"      "HILPDA"
 [806] "IGHM"        "CD48"        "FHIT"        "CD2"         "WNT9A"
 [811] "SELL"        "COL4A1"      "FCRLA"       "IL2RB"       "ITGA3"
 [816] "IL27RA"      "DLL4"        "COL16A1"     "CD52"        "CD27"
 [821] "VCAN"        "CLEC5A"      "VEGFB"       "P2RX5"       "MRC1"
 [826] "KRT10"       "AZU1"        "MGP"         "STMN1"       "DUSP2"
 [831] "CD22"        "MS4A1"       "CAMP"        "PTGIS"       "SLA"
 [836] "ACE"         "FYB1"        "TIE1"        "SLC2A1"      "TLR8"
 [841] "JAK2"        "S100B"       "XKR4"        "TNFRSF10D"   "TLR3"
 [846] "MMP9"        "COL8A1"      "LEP"         "CCL26"       "PRTN3"
 [851] "CLCF1"       "CDKN3"       "NOD2"        "CLEC12A"     "MSMB"
 [856] "CDH5"        "IL1RAP"      "RAMP2"       "NOTCH2"      "AREG"
 [861] "UPK3A"       "FZD4"        "EPHA7"       "LAMA4"       "TGFB2"
 [866] "FLT3LG"      "MKI67"       "CXCR4"       "NELL2"       "IL17RE"
 [871] "WNT7A"       "LUM"         "JAG1"        "TACSTD2"     "FGF7"
 [876] "REG1A"       "IL11RA"      "CXCR1"       "KRT13"       "FABP4"
 [881] "TCL1A"       "WNT11"       "HMGB2"       "NLRP3"       "NUSAP1"
 [886] "EMP3"        "GZMH"        "STAT4"       "CD3E"        "TPSAB1"
 [891] "CTSW"        "ADIPOQ"      "CPB1"        "DDR1"        "CD300A"
 [896] "JCHAIN"      "TWIST2"      "FAP"         "CPA3"        "VWF"
 [901] "TNNT2"       "TFEB"        "RGS5"        "IL12B"       "MZB1"
 [906] "KLRB1"       "ADGRG1"      "CCRL2"       "TGFB1"       "KRT23"
 [911] "ARID5B"      "TNFRSF11A"   "PARP1"       "TAP2"        "MIR4435-2HG"
 [916] "GPER1"       "CX3CL1"      "ANKRD1"      "MECOM"       "ITM2A"
 [921] "CEACAM6"     "CD53"        "IGHG2"       "GZMB"        "BEST1"
 [926] "CLEC7A"      "NLRC5"       "HTT"         "CSK"         "CTSG"
 [931] "IL2RA"       "IL2"         "CENPF"       "CRIP1"       "ITGB8"
 [936] "NOSIP"       "GCG"         "NRXN1"       "IFNA1"       "CCL8"
 [941] "ITGA9"       "NPR1"        "VPREB3"      "ALOX5AP"     "NRXN3"
 [946] "TLR1"        "ENTPD1"      "BMPR1A"      "CD83"        "FGF18"
 [951] "CD69"        "FGFR1"       "TNFSF15"     "TAP1"        "ADGRV1"
 [956] "IFNG"        "CDH11"       "EPCAM"       "GPR183"      "TLR7"
 [961] "CD36"        "SST"         "INHBB"       "CSF2"        "COL12A1"
 [966] "PTGS1"       "IL34"        "PGF"         "MMP14"       "CCL18"
 [971] "PTGS2"       "BMP3"        "COL1A2"      "ANGPT1"      "EOMES"
 [976] "ADGRF3"      "MYL7"        "ITGA2"       "RAG1"        "IL36G"
 [981] "NTRK2"       "VEGFC"       "MYH11"       "TGFBR1"      "KRT5"
 [986] "GZMA"        "TNFSF8"      "TBX21"       "OAS3"        "RGS13"
 [991] "TNFRSF17"    "HAVCR2"      "TNFRSF19"    "CLEC14A"     "TWIST1"
 [996] "MB"          "IDO1"        "MYH6"        "BMP7"        "ITGB6"
>
>


"""

# get the result for PCA
pca <- tiledb_scdataset$somas$RNA$obsm$members$dimreduction_pca$to_matrix()
pca[1:4,1:4]

neighbs <- tiledb_scdataset$somas$RNA$obsp$members$graph_RNA_pca_nn$to_seurat_graph()
neighbs[1:4,1:25]

## Read data into Seurat 

# Seurat object with all SOMAs 
# seurat <- tiledb_scdataset$to_seurat(batch_mode = TRUE)
# seurat
RNA_seurat <- tiledb_scdataset$to_seurat(somas = c("RNA"), batch_mode = TRUE)
RNA_seurat

RNA_seurat@meta.data <- metadata
Idents(RNA_seurat) <- RNA_seurat$cellType

markers <- FindMarkers(RNA_seurat, ident.1 = unique(Idents(RNA_seurat))[1],
                       logfc.threshold = 0.25, test.use = "roc",
                       only.pos = TRUE)
pdf("Vlnplot_H2AZ.pdf",width = 12, height = 8)
VlnPlot(RNA_seurat, features = "H2AZ1",
        log = TRUE, pt.size = 0.3)
dev.off()

# Spatial Mapping
# In the Seurat object, you might also have spatial coordinates that allow you to visualize gene expression spatially. 
# Your current approach focuses on differential expression analysis.

seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
seurat_obj <- AddMetaData(seurat_obj, metadata = transcriptCoords, col.name = "transcriptCoords")
pdf("spatial_H2AZ.pdf",width = 12, height = 8)
SpatialPlot(seurat_obj, features = "H2AZ1")
dev.off()


copy_seurat<-RNA_seurat
# check id to see seruat and transcriptCoords are match:

# Find common IDs
common_ids <- intersect(rownames(copy_seurat@meta.data), transcriptCoords$cell_id)
missing_ids_in_seurat <- setdiff(transcriptCoords$cell_id, rownames(copy_seurat@meta.data))
missing_ids_in_transcriptCoords <- setdiff(rownames(copy_seurat@meta.data), transcriptCoords$cell_id)

cat("Missing IDs in Seurat metadata:", missing_ids_in_seurat, "\n")
cat("Missing IDs in transcriptCoords:", missing_ids_in_transcriptCoords, "\n")

# add transcriptCoords to seruat
# Ensure transcriptCoords is ordered to match Seurat metadata
transcriptCoords_subset <- transcriptCoords[transcriptCoords$cell_id %in% rownames(copy_seurat@meta.data), ]
transcriptCoords_subset <- transcriptCoords_subset %>%
  dplyr::select(cell_id, slideID, fov, x_FOV_px, y_FOV_px, z_FOV_slice, target, CellComp) %>%
  dplyr::column_to_rownames(var = "cell_id")

# Reorder to match Seurat metadata row names
transcriptCoords_subset <- transcriptCoords_subset[rownames(copy_seurat@meta.data), ]


# Add metadata
copy_seurat <- AddMetaData(copy_seurat, metadata = transcriptCoords_subset)
head(copy_seurat@meta.data)

pca <- as.data.frame(tiledb_scdataset$somas$RNA$obsm$members$dimreduction_pca$to_matrix())
colorColumn <- "cellType"
pca$colorBy <- tiledb_scdataset$somas$RNA$obs$to_dataframe(attrs = colorColumn)[[1]]

pdf("PCA.pdf",width = 12, height = 8)
PCAPlot(copy_seurat, dims = c(1, 2))
dev.off()



copy_seurat <- AddMetaData(copy_seurat, metadata = transcriptCoords, col.name = "transcriptCoords")
pdf("spatial_H2AZ.pdf",width = 12, height = 8)
SpatialPlot(copy_seurat, features = "H2AZ1")
dev.off()













###
# Load required libraries
library(tiledb)
library(ggplot2)

# Extract cell coordinates
cell_coords <- as.data.frame(tiledb_scdataset$somas$RNA$members$X)

# Extract cell type annotations
cell_annotations <- as.data.frame(tiledb_scdataset$somas$RNA$members$obs)

# Combine the data into one dataframe
plot_data <- cbind(cell_coords, cell_annotations)

# Inspect the combined data
head(plot_data)


# Define a custom plot function
customImageDimPlot <- function(data, fov, x_col = "X1", y_col = "X2", color_col = "cell_type", title = "Cell Type Visualization", size = 0.3, alpha = 0.5, cols = "glasbey") {
  ggplot(data, aes_string(x = x_col, y = y_col, color = color_col)) +
    geom_point(size = size, alpha = alpha) +
    labs(title = title, x = "X coordinate", y = "Y coordinate") +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_color_manual(values = scales::hue_pal()(length(unique(data[[color_col]]))))
}

# Assuming 'cell_type' is the column for cell type annotations in plot_data
# Generate the plot
cell_type_plot <- customImageDimPlot(plot_data, fov = "lung5.rep1")

# Print the plot
print(cell_type_plot)

# Save the plot to a PDF if needed
ggsave("Cell_Type_Visualization.pdf", plot = cell_type_plot, width = 12, height = 8)

###########
###########
# Specify the slide and field of view
slide <- 2
fov <- 9

# Ensure that column names match your dataset structure
# Assuming `slide_ID_numeric` and `fov` are column names in your dataset
slideName <- unique(plot_data$Run_Tissue_name[plot_data$slide_ID_numeric == slide])

# Filter coordinates for the specified slide and FOV
fovCoords <- plot_data[plot_data$slide_ID_numeric == slide & plot_data$fov == fov, ]

# Inspect the filtered data
head(fovCoords)

# Define a custom plot function
customImageDimPlot <- function(data, x_col = "X1", y_col = "X2", color_col = "cell_type", title = "Cell Type Visualization", size = 0.3, alpha = 0.5, cols = "glasbey") {
  ggplot(data, aes_string(x = x_col, y = y_col, color = color_col)) +
    geom_point(size = size, alpha = alpha) +
    labs(title = title, x = "X coordinate", y = "Y coordinate") +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_color_manual(values = scales::hue_pal()(length(unique(data[[color_col]]))))
}

# Assuming 'cell_type' is the column for cell type annotations in fovCoords
# Generate the plot
cell_type_plot <- customImageDimPlot(fovCoords)

# Print the plot
print(cell_type_plot)

# Save the plot to a PDF if needed
ggsave("Cell_Type_Visualization.pdf", plot = cell_type_plot, width = 12, height = 8)



