#' Interchangeable single-cell data formats
#'
#' Convert a Seurat object to a scanpy object.
#'
#' This R function accepts a Seurat object and write a scanpy h5ad file to the disk.
#'
#' multiple parameters could be specified to control the conversion behaviour.
#' @param seu \code{seu} is the seurat object to be converted
#' @param save_name  \code{save_name} is the prefix of the output file's name.
#' @param save_embedding \code{save_embedding} shoule be set to TRUE when you
#' hope put embeddings (e.g. PCA, UMAP) into the converted file.
#' @return No value returned. an h5ad file is written to the disk.
#' @author Sijie Chen
#' @export
CreateScanpyObject<-function(seu, save_name, save_embedding = TRUE){
    sample <- save_name
    if (save_embedding){
        tictoc::tic("saving embeddings")
        dr.umap <- Seurat::Embeddings(object=seu,reduction ="umap")[rownames(seu@meta.data), 1:2]
        dr.pca  <- Seurat::Embeddings(object=seu,reduction ="pca") [rownames(seu@meta.data) ,1:2]
        seu@meta.data <- cbind(seu@meta.data,dr.umap, dr.pca)        
        tictoc::toc()
    }


    tictoc::tic("extracting dge,var,obs")
    dge <- Seurat::GetAssayData(seu, assay = "RNA", slot = "counts") # gene by cell
    var <- rownames(dge) # rownames are genes
    obs <- seu@meta.data
    tictoc::toc()


    tictoc::tic("writing dge,var,obs")
    tempFolder = tempdir(check = TRUE)
    dgePath = file.path(tempFolder, paste(sample, "_dge.spmtx" ,sep=""))
    varPath = file.path(tempFolder, paste(sample, "_var.tsv" ,sep=""))
    obsPath = file.path(tempFolder, paste(sample, "_obs.tsv" ,sep=""))
    
    #library(Matrix)
    Matrix::writeMM(dge, file = dgePath)
    write.table(var, sep = "\t" ,row.names =F, col.names = F, file=varPath)
    write.table(obs, sep = "\t", row.names =T, col.names = T, file=obsPath)
    tictoc::toc()
    
    print(dgePath)
    print(varPath)
    print(obsPath)
    
    tictoc::tic("converting to scanpy")
    build_anndata_script=system.file("python/build_anndata.py",package="scinterchange")
    message(paste("use python=",build_anndata_script))
    reticulate::source_python(file = build_anndata_script)
    build_anndata_from_files(dgePath=dgePath,
                             varPath=varPath,
                             obsPath=obsPath,
                             output_file=paste(sample, ".h5ad" ,sep=""),
                             with_embedding=save_embedding)
    tictoc::toc()
}