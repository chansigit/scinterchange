def build_anndata_from_files(dgePath, varPath, obsPath, output_file, with_embedding=True):
    import pandas as pd
    from scipy.io import mmread
    sparse = mmread(dgePath)
    X = sparse.T.tocsc()

    var = pd.DataFrame(index=pd.read_csv(varPath,delimiter='\t',
                                         encoding='utf-8', header=None,
                                         names=["gene"]).to_numpy().flatten() )
    obs = pd.read_csv(obsPath, delimiter='\t', encoding='utf-8')

    import anndata as ad
    adata = ad.AnnData(X, obs=obs, var=var, dtype='int32')
    adata.var_names_make_unique()
    if with_embedding:
        import numpy as np
        adata.obsm["X_umap"]=np.array(adata.obs[["UMAP_1","UMAP_2"]])
        adata.obsm["X_pca"]=np.array(adata.obs[["PC_1","PC_2"]])
    adata.write_h5ad(filename=output_file)