# --- Spatial Transcriptomics Analysis Pipeline ---
# Description:
# This script performs spatial transcriptomics analysis on CHC23 tissue data (GSE245908).
# It includes preprocessing, quality control, clustering, and visualization.

# Required Libraries
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
import os
import seaborn as sns
import gzip
import shutil

# --- File Paths ---
# Input compressed file (tissue positions)
input_path = '/home/jiew/GSE245908/GSM7850823_CHC23_tissue_positions_list.csv.gz'

# Output decompressed file
output_path = '/home/jiew/GSE245908/GSM7850823_CHC23_tissue_positions_list.csv'

# Decompress the input file
with gzip.open(input_path, 'rb') as f_in:
    with open(output_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

# --- Set Working Directory ---
# Change to results directory
os.chdir('/home/jiew/GSE245908/CHC23_result')

# --- Load and Inspect Data ---
# Load Visium spatial transcriptomics data
adata1 = sc.read_visium("/home/jiew/GSE245908/CHC23_Processed")
adata1.var_names_make_unique() # Ensure unique gene names
adata1 # Print summary of the AnnData object
adata1.obsm['spatial'].shape # Print summary of the AnnData object

# --- Visualize Raw Spatial Data ---
adata1.obs['thing'] = 'a' # Dummy variable for initial visualization
plt.rcParams['figure.figsize'] = (8,8)
sc.pl.spatial(adata1,color ='thing')
plt.savefig('check.png')
plt.close()  

# --- Quality Control ---
# Identify mitochondrial genes
adata1.var["mt-"] = adata1.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata1, qc_vars=["mt-"], inplace = True)
print(adata1.obs.head())  # View quality metrics

# Plot quality control metrics
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.distplot(adata1.obs["total_counts"], kde=False, ax=axs[0])
sns.distplot(adata1.obs["total_counts"][adata1.obs["total_counts"] <10000], kde=False, bins =40, ax=axs[1])
sns.distplot(adata1.obs["n_genes_by_counts"], kde=False, bins =60,ax=axs[2])
sns.distplot(adata1.obs["n_genes_by_counts"][adata1.obs["n_genes_by_counts"] <4000], kde=False, bins =40, ax=axs[3])
plt.show()
plt.savefig('qc_metrics.png')
plt.close()  


# --- Filter Cells and Genes ---
# Filter cells based on total counts
sc.pp.filter_cells(adata1, min_counts =1000)
sc.pp.filter_cells(adata1, max_counts =80000)
sc.pl.violin(adata1, ["pct_counts_mt-"], jitter = 0.4)
plt.savefig('qc.png')
plt.close()  

# Filter cells with high mitochondrial gene expression
adata1 = adata1[adata1.obs["pct_counts_mt-"] < 28]

# Filter genes based on cell expression
sc.pp.filter_genes(adata1, min_cells =3)

plt.rcParams['figure.figsize'] = (8,8)
sc.pl.spatial(adata1,color ='thing')
plt.savefig('check_filtered_a.png')
plt.close()  


# --- Normalization and Log Transformation ---
sc.pp.normalize_total(adata1,inplace=True) # Normalize data
sc.pp.log1p(adata1) # Log-transform data

# --- Identify Highly Variable Genes ---
sc.pp.highly_variable_genes(adata1,flavor="seurat", n_top_genes=2000)

# --- Dimensionality Reduction and Clustering ---
sc.pp.pca(adata1)
sc.pp.neighbors(adata1)
sc.tl.umap(adata1)
sc.tl.leiden(adata1)

# --- Visualize UMAP Clusters ---
plt.rcParams["figure.figsize"] = (4,4)
sc.pl.umap(adata1, color =["total_counts", "n_genes_by_counts","leiden"], wspace =0.4)
plt.savefig('leiden_a.png')
plt.close()  

plt.rcParams["figure.figsize"] = (4,4)
sc.pl.spatial(adata1, img_key="hires",color =["total_counts", "n_genes_by_counts"])
plt.savefig('spatial_a_.png')
plt.close() 

sc.pl.spatial(adata1, img_key="hires",color ="leiden", size=1.5)
plt.savefig('spatial_leiden_a_a.png')
plt.close() 

# --- Differential Expression Analysis ---
# Rank genes for each cluster using Wilcoxon rank-sum test
sc.tl.rank_genes_groups(adata1,"leiden",method='wilcoxon')

# Visualize top marker genes in heatmap
sc.pl.rank_genes_groups_heatmap(adata1, n_genes=10)
plt.savefig('heatmap.png')
plt.close() 

# Extract ranked genes into a DataFrame
results = adata1.uns["rank_genes_groups"]
('0','1','2','3','4')
out =np.array([[0,0,0,0,0]])
for group in results['names'].dtype.names:
    out = np.vstack((out, np.vstack((results['names'][group],
                                     results['scores'][group],
                                     results['pvals_adj'][group],
                                     results['logfoldchanges'][group],
                                     np.array([group]*len(results['names'][group])).astype('object'))).T))

markers = pd.DataFrame(out[1:], columns=['Gene','scores','pval_adj','lfc','cluster'])

# Filter significant markers
markers = markers[(markers.pval_adj <0.05)&abs(markers.lfc>1)]
# markers[markers.cluster =='2']
markers.to_csv('markers.csv', index=False)


# --- Visualize Marker Genes ---
sc.pl.spatial(adata1, img_key='hires', color =['ACTA2','TIMP1','COL1A1'],size =1.5,)
plt.savefig('Fibrosis_Cell_.png')
plt.close() 
quit()

# --- Gene Expression Dotplot ---
genes_of_interest = ['PTPRC','CD3E','CD3D','KLRB1','CD4','CD8A', 'ACTA2', 'TIMP1', 'COL1A1']
sc.pl.dotplot(adata1, genes_of_interest, groupby='leiden', dot_max=0.5, standard_scale='var')
plt.savefig('Bubble_plot_a.png')
plt.close() 

 

# --- Re-clustering with Different Resolution ---
sc.tl.leiden(adata1, resolution=1.0)
# Save the cluster labels to adata
adata1.obs['leiden'] = adata1.obs['leiden'].astype('category')
adata1

# --- Save Processed Data ---
adata1.write('CHC23_adata.h5ad')