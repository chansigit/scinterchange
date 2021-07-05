# scinterchange

## introduction

There are many single-cell data formats. We usually need to convert them from one format to another. 

**scinterchange** supports conversions:

- from R seurat to python scanpy
- ...
- 

contact: Sijie Chen (chansigit@gmail.com )

---

## installation

- in China `devtools::install_git("https://gitee.com/chansigit/scinterchange.git")`

- global `devtools::install_git("https://github.com/chansigit/scinterchange.git")`

---

## Use

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

