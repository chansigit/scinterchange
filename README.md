# scinterchange

## introduction

There are many single-cell data formats. We usually need to convert them from one format to another. 

**scinterchange** supports conversions:

- from R seurat to python scanpy
- from python scanpy to R seurat
- 

contact: Sijie Chen (chansigit@gmail.com )

---

## installation
- Now I do not recommended my scinterchange package, please see the pure reticulate method below (credit to Theis Lab)

- in China `devtools::install_git("https://gitee.com/chansigit/scinterchange.git")`

- global `devtools::install_git("https://github.com/chansigit/scinterchange.git")`

---

## Use


### Convert a Seurat object to a scanpy object (pure reticulate, recommended)
```R
# solution from https://theislab.github.io/scanpy-in-R
library(reticulate)
sc <- import("scanpy")
use_python("/home/hcauser/anaconda3/envs/r410py37/bin/python", required = T)

adata_seurat <- sc$AnnData(
    X   = t(GetAssayData(seurat)),
    obs = seurat[[]],
    var = GetAssay(seurat)[[]]
)
adata_seurat$obsm$update(umap = Embeddings(seurat, "umap"))
sp <- import("scipy.sparse")
adata_seurat$X  <- sp$csc_matrix(adata_seurat$X)
adata_seurat

adata_seurat$write_h5ad("*****.h5ad")
```

### Convert a scanpy object to a Seurat object (pure reticulate, recommended)
```R
# solution from https://theislab.github.io/scanpy-in-R
# first save scanpy object in the h5ad format in python environment
library(reticulate)
sc <- import("scanpy")
adata <- sc$read_h5ad("*******.h5ad")
adata

## NOTE that if adata$X is a sparse matrix, t() will report an error. 
## consider densify the matrix first and then transpose
## `t(as.matrix(adata$X))`
## or use t_sp() provided by R package bayNorm if the matrix is really big
## # BiocManager::install("bayNorm")
## exprs <- bayNorm::t_sp(adata$X)
exprs <- t(adata$X)
colnames(exprs) <- adata$obs_names$to_list()
rownames(exprs) <- adata$var_names$to_list()

# Create the Seurat object
seurat <- CreateSeuratObject(exprs)

# Set the expression assay
seurat <- SetAssayData(seurat, "data", exprs)

# Add observation metadata
seurat <- AddMetaData(seurat, adata$obs)

# Add fetaure metadata
seurat[["RNA"]][["n_cells"]] <- adata$var["n_cells"]

# Add embedding
embedding <- adata$obsm["X_umap"]
rownames(embedding) <- adata$obs_names$to_list()
colnames(embedding) <- c("umap_1", "umap_2")
seurat[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")
```

### Convert a Seurat object to a scanpy object

```R
library(reticulate)
use_python("/home/hcauser/anaconda3/envs/r410py37/bin/python", required = T)

CreateScanpyObject(seu, save_name="ilcTemp", save_embedding=TRUE)

#saving embeddings: 0.007 sec elapsed
#extracting dge,var,obs: 0.014 sec elapsed
#writing dge,var,obs: 1.042 sec elapsed
#[1] "/tmp/RtmppZhWii/ilcTemp_dge.spmtx"
#[1] "/tmp/RtmppZhWii/ilcTemp_var.tsv"
#[1] "/tmp/RtmppZhWii/ilcTemp_obs.tsv"
#converting to scanpy: 3.623 sec elapsed
# 
#you will see an ouput file as ilcTemp.h5ad in your working directory.
```
