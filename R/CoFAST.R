
get.indexList <- function(seu, batch) {
  if (!batch %in% colnames(seu@meta.data)) {
    stop("batch is not in the columns of meta.data!")
  }
  batch_ID <- unique(seu@meta.data[, batch])
  indexList <- lapply(batch_ID, function(ID) {
    which(seu@meta.data[, batch] == ID)
  })
  return(indexList)
}


compute.AdjList <- function(
  seu, spatialCoords, batch = NULL,
  type = c("fixed_distance", "fixed_number"),
  platform = c("Others", "ST", "Visium"), radius_vec = NULL, ...) {
  calAdjd <- function(coord.mat, k = 6) {
    n <- nrow(coord.mat)
    nn.idx <- RANN::nn2(data = coord.mat, k = k + 1)$nn.idx
    j <- rep(1:n, each = k + 1)
    i <- as.vector(t(nn.idx))
    Adj <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(n, n))
    diag(Adj) <- 0
    return(Adj)
  }

  if (is.null(batch)) {
    idxList <- list(seq_len(ncol(seu)))
  } else {
    idxList <- get.indexList(seu, batch)
  }
  if (length(radius_vec) == 1) {
    rep(radius_vec, length(idxList))
  }
  message("Calculate the adjacency matrix...")
  posList <- lapply(idxList, function(idx) {
    as.matrix(seu@meta.data[idx, spatialCoords])
  })
  type <- match.arg(type)
  platform <- match.arg(platform)
  d <- length(spatialCoords)
  if (d == 2) {
    if (tolower(type) == "fixed_distance") {
      if (tolower(platform) %in% c("st", "visium")) {
        AdjList <- pbapply::pblapply(
          posList, PRECAST::getAdj_reg, platform = platform)
      } else {
        if (is.null(radius_vec)) {
          AdjList <- pbapply::pblapply(posList, function(x) {
            DR.SC::getAdj_auto(x, ...)
          })
        } else {
          AdjList <- pbapply::pblapply(
            seq_along(posList), function(r) {
              DR.SC::getAdj_auto(
                posList[[r]], radius.upper = radius_vec[r])
            })
        }
      }
    } else if (tolower(type) == "fixed_number") {
      AdjList <- pbapply::pblapply(posList, PRECAST::getAdj_fixedNumber, ...)
    } else {
      stop("AddAdjList: Unsupported adjacency  type \"", type, "\".")
    }
  } else {
    AdjList <- pbapply::pblapply(posList, calAdjd, ...)
  }

  seu@misc[["AdjList"]] <- AdjList
  return(seu)
}


#' Determine the dimension of low dimensional embedding
#'
#' @useDynLib ProFAST, .registration = TRUE
#' @description
#' This function estimate the dimension of low dimensional embedding for a given cell by gene expression matrix. For more details, see Franklin et al. (1995) and Crawford et al. (2010).
#' @references
#' 1. Franklin, S. B., Gibson, D. J., Robertson, P. A., Pohlmann, J. T., & Fralish, J. S. (1995). Parallel analysis: a method for determining significant principal components. Journal of Vegetation Science, 6(1), 99-106.
#' 
#' 2. Crawford, A. V., Green, S. B., Levy, R., Lo, W. J., Scott, L., Svetina, D., & Thompson, M. S. (2010). Evaluation of parallel analysis methods for determining the number of factors.Educational and Psychological Measurement, 70(6), 885-901.
#'
#' @param object A Seurat or matrix object
#' @param ... Arguments passed to other methods
#' @rdname diagnostic.cor.eigs
#' @export diagnostic.cor.eigs
#'
#' @return A data.frame with attribute `q_est` and `plot`, which is the estimated dimension of low dimensional embedding. In addition, this data.frame containing the following components:
#' \itemize{
#'   \item q - The index of eigen values.
#'   \item eig_value - The eigen values on observed data.
#'   \item eig_sim - The mean value of eigen values of n.sims simulated data.
#'   \item q_est - The selected dimension in attr(obj, 'q_est').
#'   \item plot - The plot saved in attr(obj, 'plot').
#' }
#'
diagnostic.cor.eigs <- function(object, ...) {  
  UseMethod("diagnostic.cor.eigs", object = object)  
}

#' Determine the dimension of low dimensional embedding
#' @importFrom irlba irlba
#' @importFrom furrr future_map
#' @importFrom future plan
#' @importFrom stats rnorm 
#' @import ggplot2
#'
#' @param q_max the upper bound of low dimensional embedding. Default is 50.
#' @param plot a indicator of whether plot eigen values.
#' @param n.sims number of simulaton times. Default is 10.
#' @param parallel a indicator of whether use parallel analysis.
#' @param ncores the number of cores used in parallel analysis. Default is 10.
#' @param seed a postive integer, specify the random seed for reproducibility
#'
#'
#' @rdname diagnostic.cor.eigs
#'
#' @export
#'
#' @examples
#' n <- 100
#' p <- 50
#' d <- 15
#' object <- matrix(rnorm(n*d), n, d) %*% matrix(rnorm(d*p), d, p)
#' diagnostic.cor.eigs(object, n.sims=2)
diagnostic.cor.eigs.default <- function(
  object,  q_max = 50,  plot = TRUE,
  n.sims = 10, parallel = TRUE, ncores = 10, seed=1, ...) {
  if (!is.numeric(q_max)) {
    stop("q_max should be a numeric!")
  }
  if (q_max <= 0) {
    stop("q_max should larger than 0!")
  }
  if (q_max >= ncol(object)) {
    warning("q_max is not less than the number of columns of X. Set it as ncol(X) - 1")
    q_max <- ncol(object) - 1
  }
  Y <- scale(object)
  n <- nrow(Y)
  p <- ncol(Y)
  svdX <- irlba::irlba(A = Y / sqrt(n), nv = q_max)

  dvec <- svdX$d^2 ## is the eigenvalues of correlation matrix
  # if ((!dir.exists(dir_name)) && (any(c(plot, save_eigen)))) {
  #   dir.create(dir_name)
  # }
  # if (save_eigen) {
  #   save(dvec, file = paste0("./", dir_name, '/cor_eigs_values.rds'))
  # }

  ### simulate data
  corr_fun <- function(i, n, p) {
    #set.seed(i)
    X1 <- matrix(rnorm(n * p), n, p)
    svdX1 <- irlba::irlba(A = X1 / sqrt(n), nv = q_max)
    return(svdX1$d^2)
  }
  if (parallel) {
    #library(furrr)
    #library(future)
    future::plan('multicore', workers = ncores)
    eig.mat.sim <- furrr::future_map(1:n.sims,corr_fun, n = n, p = p, .progress = TRUE, .options = furrr::furrr_options(seed = seed))
    eig.mat.sim <- Reduce(cbind, eig.mat.sim)
  } else {
    eig.mat.sim <- pbapply::pbsapply(1:n.sims, corr_fun, n = n, p = p)
  }
  dat <- data.frame(
    'q' = seq_along(dvec),
    'eig_value' = dvec,
    'eig_sim' = rowMeans(eig.mat.sim))
  q_s <- which(dat$eig_sim > dat$eig_value)[1] - 1
  if (plot) {
    p1 <- ggplot(data = dat, aes_string(x = 'q', y = 'eig_value')) +
      geom_line(linewidth = 1) +
      geom_point(size = 1.5) +
      geom_line(aes_string(x = 'q', y = 'eig_sim'), color = "red")+
      geom_vline(xintercept = q_s) +
      theme_classic(base_size = 16)
    print(p1)
    attr(dat, "plot") <- p1
    # write_fig(
    #   p1, filename = "parallel_analysis_plot", dir_name = dir_name)
  }
  if (is.na(q_s)) {
    q_s <- q_max
  }
  attr(dat, "q_est") <- q_s
  return(dat)
}






#' @param assay an optional string, specify the name of assay in the Seurat object to be used.
#' @param slot an optional string, specify the name of slot.
#' @param nfeatures an optional integer, specify the number of features to select as top variable features. Default is 2000.
#' @param ... Other arguments passed to \code{\link{diagnostic.cor.eigs.default}}.
#'
#' @method diagnostic.cor.eigs Seurat
#' @rdname diagnostic.cor.eigs
#' @importFrom Seurat DefaultAssay
#' @importFrom Matrix t
#'
#' @export
#'
diagnostic.cor.eigs.Seurat <- function(
  object, assay = NULL, slot = "data", nfeatures = 2000, q_max = 50,
  seed = 1,...) {
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  X_all <- Seurat::GetAssayData(object = object, slot = slot, assay = assay)
  if (length(object@assays[[assay]]@var.features) == 0) {
    object <- Seurat::FindVariableFeatures(
      object, nfeatures = min(nfeatures, nrow(object)), assay = assay)
  }
  var.features <- object@assays[[assay]]@var.features
  X <- Matrix::t(X_all[var.features, ])
  rm(X_all)
  n <- ncol(X)
  if (n > 10000) {
    set.seed(seed)
    X <- X[sample(1:n, 10000), ]
  }
  keep_idx <- (Matrix::colSums(X) != 0)
  X <- X[, keep_idx]
  dat_cor <- diagnostic.cor.eigs(X, q_max = q_max, ...)
  dat_cor
}

#' Cell-feature coembedding for scRNA-seq data
#' @description Cell-feature coembedding for scRNA-seq data based on FAST model.
#' @param object a Seurat object.
#' @param assay an optional string, specify the name of assay in the Seurat object to be used, `NULL` means default assay in seu.
#' @param slot an optional string, specify the name of slot.
#' @param nfeatures an optional integer, specify the number of features to select as top variable features. Default is 2000.
#' @param q an optional positive integer, specify the dimension of low dimensional embeddings to compute and store. Default is 10.
#' @param reduction.name an optional string, specify the dimensional reduction name, `ncfm` by default.
#' @param weighted an optional logical value, specify whether use weighted method.
#' @param var.features an optional string vector, specify the variable features used to calculate cell embedding.
#'
#' @rdname NCFM
#'
#' @importFrom Seurat DefaultAssay GetAssayData CreateDimReducObject FindVariableFeatures
#'
#' @export
#'
#' @examples
#' data(pbmc3k_subset)
#' pbmc3k_subset <- NCFM(pbmc3k_subset)
NCFM <- function(
  object, assay = NULL, slot = "data", nfeatures = 2000, q = 10,
  reduction.name = "ncfm", weighted = FALSE, var.features = NULL) {
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  X_all <- as.matrix(Seurat::GetAssayData(
    object = object, slot = slot, assay = assay))
  if (is.null(var.features)) {
    if (length(object@assays[[assay]]@var.features) == 0) {
      stop("NCFM: please find the variable features using Seurat::FindVariableFeatures or DR.SC::FindSVGs before running this function!")
    }
    var.features <- object@assays[[assay]]@var.features
  } else {
    var.features <- intersect(
      var.features, rownames(X_all)) # slot is scale.data, restrict the var.features to be the features in scale.data
  }
  if(slot != 'data'){
    X_data <- as.matrix(Seurat::GetAssayData(
      object = object, slot = 'data', assay = assay))
  }else{
    X_data <- X_all
  }
  res <- Factor_nc(
    X_slot = X_all[var.features,], X_data=X_data,  q = q, reduction.name = reduction.name, weighted = weighted,
    features = NULL)
  cellsCoordinates <- res$cellsCoordinates
  featuresCoordinates <- res$featuresCoordinates

  object@reductions[[reduction.name]] <- Seurat::CreateDimReducObject(
    embeddings = cellsCoordinates,
    loadings = featuresCoordinates[rownames(object), ],
    key = paste0(reduction.name, "_"), assay = assay)
  return(object)
}



#' @importFrom Matrix t
#' @importFrom irlba irlba
Factor_nc <- function(
  X_slot, X_data, q = 10, reduction.name = "Fac", weighted = FALSE, features = NULL) {
  if (q <= 1) {
    stop("q must be greater than or equal to 2!")
  }
  if (is.null(features)) {
    features <- row.names(X_slot)
  } else {
    features <- intersect(features, row.names(X_slot))
  }
  X_slot <- as.matrix(X_slot)
  tstart <- Sys.time()
  px <- length(features)
  if (weighted) {
    fit <- irlba::irlba(
      A = Matrix::t(X_slot[features, ]), nu = q, nv = 1, work = sqrt(px * q))
    ce_cell <- fit$u %*% diag(fit$d[1:q])
  } else {
    ce_cell <- irlba::irlba(
      A = Matrix::t(X_slot[features, ]), nu = q, nv = 1, work = sqrt(px * q))$u
  }

  ce_gene <- gene_embed_cpp(X_data, ce_cell)

  component <- paste0(reduction.name, "_", seq_len(ncol(ce_cell)))
  colnames(ce_cell) <- component
  row.names(ce_cell) <- colnames(X_data)
  colnames(ce_gene) <- component
  row.names(ce_gene) <- row.names(X_data)
  .logDiffTime(sprintf(paste0("%s Finish CoFAST"), "*****"), t1 = tstart, verbose=TRUE)
  output <- list(
    cellsCoordinates = ce_cell,
    featuresCoordinates = ce_gene
  )
  return(output)
}


#' Cell-feature coembedding for SRT data
#' @description Run cell-feature coembedding for SRT data based on FAST model.
#' @param object a Seurat object.
#' @param Adj_sp a sparse matrix, specify the adjacency matrix among spots.
#' @param assay an optional string, the name of assay used.
#' @param slot an optional string, the name of slot used.
#' @param nfeatures an optional postive integer, the number of features to select as top variable features. Default is 2000.
#' @param q an optional positive integer, specify the dimension of low dimensional embeddings to compute and store. Default is 10.
#' @param reduction.name an optional string, dimensional reduction name, `fast` by default.
#' @param var.features an optional string vector, specify the variable features, used to calculate cell embedding.
#' @param ... Other argument passed  to the \code{\link{FAST_run}}.
#'
#' @importFrom Seurat DefaultAssay GetAssayData CreateDimReducObject
#'
#' @export
#'
#' @examples
#' data(CosMx_subset)
#' pos <- as.matrix(CosMx_subset@meta.data[,c("x", "y")])
#' Adj_sp <- AddAdj(pos)
#' # Here, we set maxIter = 3 for fast computation and demonstration.
#' CosMx_subset <- NCFM_fast(CosMx_subset, Adj_sp = Adj_sp, maxIter=3)
#' 
NCFM_fast <- function(
  object, Adj_sp, assay = NULL, slot = "data", nfeatures = 2000, q = 10,
  reduction.name = "fast", var.features = NULL, ...) {
  
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  tstart <- Sys.time()
  X_all <- (Seurat::GetAssayData(
    object = object, slot = slot, assay = assay)) # as.matrix
  if (is.null(var.features)) {
    if (length(object@assays[[assay]]@var.features) == 0) {
      stop("NCFM: please find the variable features using Seurat::FindVariableFeatures or DR.SC::FindSVGs before running this function!")
    }
    var.features <- object@assays[[assay]]@var.features
  } else {
    var.features <- intersect(
      var.features, rownames(X_all))
  }
  if(slot != 'data'){
    X_data <- as.matrix(Seurat::GetAssayData(
      object = object, slot = 'data', assay = assay)) ## ensure X_data always includes all features in seu.
  }else{
    X_data <- X_all 
  }
  
  res <- FAST_single_coembed(
    X_slot = X_all[var.features,], X_data = X_data, Adj_sp = Adj_sp, q = q, var.features = NULL, ...)
  cellsCoordinates <- res$cellsCoordinates
  featuresCoordinates <- res$featuresCoordinates
  colnames(cellsCoordinates) <- paste0("fast", 1:q)
  .logDiffTime(sprintf(paste0("%s Finish CoFAST"), "*****"), t1 = tstart, verbose=TRUE)
  feature_names <- intersect(row.names(featuresCoordinates), rownames(object))
  object@reductions[[reduction.name]] <- Seurat::CreateDimReducObject(
    embeddings = cellsCoordinates[colnames(object), ],
    loadings = featuresCoordinates[feature_names, ],
    key = paste0(gsub("_", "", reduction.name), "_"), assay = assay)

  return(object)
}


FAST_single_coembed <- function(
  X_slot, X_data, Adj_sp, q = 15, var.features = NULL, ...) {
  if (is.null(var.features)) {
    var.features <- row.names(X_slot)
  } else {
    var.features <- intersect(var.features, row.names(X_slot))
  }

  reslist <- FAST_run(
    XList = list(Matrix::t(X_slot[var.features, ])), AdjList = list(Adj_sp),
    q = q, fit.model = "gaussian", ...)
  ce_cell <- reslist$hV[[1]]
  row.names(ce_cell) <- colnames(X_slot)
  rm(reslist)

  ce_gene <- gene_embed_cpp(as.matrix(X_data), ce_cell)

  row.names(ce_gene) <- row.names(X_data)

  output <- list(
    cellsCoordinates = ce_cell,
    featuresCoordinates = ce_gene
  )
  return(output)
}



pdistance.matrix <- function (Ar, Br, eta = 1e-10) {
  dis <- pdistance_cpp(Ar, Br, eta = eta)
  rownames(dis) <- rownames(Ar)
  colnames(dis) <- rownames(Br)
  return(dis)
}

#'Calculate the cell-feature distance matrix
#' @description  Calculate the cell-feature distance matrix based on coembeddings.
#' @param object a Seurat object.
#' @param reduction a opstional string, dimensional reduction name, `fast` by default.
#' @param assay.name a opstional string, specify the new generated assay name, `distce` by default.
#' @param eta an optional postive real, a quantity to avoid numerical errors. 1e-10 by default.
#'
#' @details This function calculate the distance matrix between cells/spots and features, and then put the distance matrix in a new generated assay. This 
#' distance matrix will be used in the siganture gene identification.
#' @importFrom Seurat Loadings Embeddings CreateAssayObject
#'
#' @export
#'
#' @examples
#' data(pbmc3k_subset)
#' pbmc3k_subset <- NCFM(pbmc3k_subset)
#' pbmc3k_subset <- pdistance(pbmc3k_subset, "ncfm")
pdistance <- function(
  object, reduction = "fast", assay.name = "distce", eta = 1e-10) {
  if(!inherits(object, "Seurat")) stop("pdistance: object must be a Seurat object!")
  f_ebd <- Seurat::Loadings(object, reduction)
  c_ebd <- Seurat::Embeddings(object, reduction)
  message("Calculate co-embedding distance...")
  pdis <- pdistance.matrix(f_ebd, c_ebd, eta)
  object@assays[[assay.name]] <- Seurat::CreateAssayObject(data = pdis)
  object@assays[[assay.name]]@misc$reduction <- reduction
  object@assays[[assay.name]]@key <- paste0(assay.name, "_")
  object
}



# Find signature genes ----------------------------------------------------

#' Obtain the top signature genes and related information
#' @description Obtain the top signature genes and related information.
#' @param df.list a list that is obtained by the function \code{\link{find.signature.genes}}.
#' @param ntop an optional positive integer, specify the how many top signature genes extracted, default as 5.
#' @param expr.prop.cutoff an optional postive real ranging from 0 to 1,  specify cutoff of expression proportion of  features, default as 0.1.
#' @return return  a `data.frame` object with four columns: `distance`,`expr.prop`, `label` and `gene`.
#' @details Using this funciton, we obtain the top signature genes and organize them into a data.frame. The `row.names` are gene names.
#' The colname `distance` means the distance between gene (i.e., VPREB3) and cells with the specific cell type (i.e., B cell),
#'  which is calculated based on the coembedding of genes and cells in the coembedding space. The distance is smaller, the association between gene and the cell type is stronger. 
#'  The colname `expr.prop` represents the expression proportion of the gene (i.e., VPREB3) within the cell type (i.e., B cell). 
#'  The colname `label` means the cell types and colname `gene` denotes the gene name. 
#'  By the data.frame object, we know `VPREB3` is the one of the top signature gene of B cell.
#' @seealso None
#' @references None
#' @export
#' @importFrom  pbapply pblapply
#' @importFrom utils head
#' @examples
#' library(Seurat)
#' data(pbmc3k_subset)
#' pbmc3k_subset <- pdistance(pbmc3k_subset, reduction='ncfm')
#' df_list_rna <- find.signature.genes(pbmc3k_subset)
#' dat.sig <- get.top.signature.dat(df_list_rna, ntop=5)
#' head(dat.sig)

get.top.signature.dat <- function(df.list, ntop=5, expr.prop.cutoff=0.1) {
  top5 <- pbapply::pblapply(seq_along(df.list), function(j) {
    #message("j = ", j)
    #j <- 1
    x <- df.list[[j]]
    df_co_sub <- x[x$expr.prop > expr.prop.cutoff, ]
    if(nrow(df_co_sub) == 0) df_co_sub <- x
    df_co_sub$label <- names(df.list)[j]
    df_co_sub$gene <- row.names(df_co_sub)
    head(df_co_sub[order(df_co_sub$distance), ], ntop)
  })
  dat <- Reduce(rbind, top5)
  return(dat)
}


#' Find the signature genes for each group of cell/spots
#' @description Find the signature genes for each group of cell/spots based on coembedding distance and expression ratio.
#' @param seu a Seurat object with coembedding in the reductions slot wiht component name reduction.
#' @param distce.assay an optional character, specify the assay name that constains distance matrix beween cells/spots and features, default as `distce` (distance of coembeddings).
#' @param ident an optional character in columns of metadata,  specify the group of cells/spots. Default as NULL, use Idents as the group.
#' @param expr.prop.cutoff an optional postive real ranging from 0 to 1,  specify cutoff of expression proportion of  features, default as 0.1.
#' @param assay an optional character,  specify the assay in seu, default as NULL, representing the default assay in seu.
#' @param genes.use an optional string vector, specify genes as the signature candidates.
#' @return return a list with each component a data.frame object having two columns: `distance` and `expr.prop`.
#' @details In each data.frame object of the returned value, the row.names are gene names, and these genes are sorted by decreasing order of `distance`. User can define the signature genes as top n genes in distance and that the `expr.prop` larger than a cutoff. We set the cutoff as 0.1.
#' @seealso None
#' @references None
#' @export
#' @importFrom  Seurat DefaultAssay Idents
#' @importFrom  pbapply pblapply
#' @examples
#' library(Seurat)
#' data(pbmc3k_subset)
#' pbmc3k_subset <- pdistance(pbmc3k_subset, reduction='ncfm')
#' df_list_rna <- find.signature.genes(pbmc3k_subset)
#'
find.signature.genes <- function(seu, distce.assay='distce', ident=NULL, expr.prop.cutoff=0.1,
                                 assay=NULL, genes.use=NULL){

  if(is.null(ident)){
    cell_label_vec <- Idents(seu)
  }else{
    if(!is.element(ident, colnames(seu@meta.data))) stop("ident must be NULL or one column name in meta.data!")
    cell_label_vec <- seu@meta.data[,ident]
  }
  if(is.null(assay)){
    assay <- DefaultAssay(seu)
  }
  cell_ID <- sort(as.character(unique(cell_label_vec)))
  barcodeList <- lapply(cell_ID, function(x){
    idx <- which(cell_label_vec== x)
    colnames(seu)[idx]
  } )
  names(barcodeList) <- cell_ID


  df_co_list <- pbapply::pblapply(seq_along(cell_ID), function(r){
    cell.set.tmp <- barcodeList[[r]]
    df_co_genes.tmp <- gene.activity.score.seu(seu, cell.set=cell.set.tmp,assay=assay,
                                               distce.assay=distce.assay, genes.use=genes.use)
    df_co_genes.tmp[df_co_genes.tmp$expr.prop>expr.prop.cutoff,]
  })
  names(df_co_list) <- cell_ID
  df_co_list <- lapply(names(df_co_list), function(x) {
    df <- df_co_list[[x]]
    df$label <- x
    df$gene <- rownames(df)
    df
  })
  names(df_co_list) <- cell_ID
  return(df_co_list)
}

#' @importFrom  Seurat DefaultAssay
#'
gene.activity.score.seu <- function(seu,  cell.set, distce.assay='distce', assay=NULL, cells.use=NULL, genes.use=NULL){

  if (is.null(genes.use)) {
    ## update by Wei Liu, log: 2024-02-24
    ## Take the features in distance matrix for signature gene candidate set.
    genes.use <- rownames(seu@assays[[distce.assay]])
  }
  if (is.null(assay)) {
    assay <- DefaultAssay(seu)
  }
  if (is.null(cells.use)) {
    cells.use <- colnames(seu)
  }
  distce <- seu@assays[[distce.assay]][genes.use, cells.use]
  if (length(cell.set) > 1) {
    genes.expr.prop <- apply(seu[[assay]]@data[genes.use, cell.set], 1, function(x) mean(x>0))
  } else {
    genes.expr.prop <- seu[[assay]]@data[genes.use, cell.set] > 0
  }
  names(genes.expr.prop) <- genes.use
  gene.activity.score(distce, genes.expr.prop, cells.use=cells.use, cell.set=cell.set)
}


gene.activity.score <- function (distce, genes.expr.prop, cells.use, cell.set){

  genes.use <- names(genes.expr.prop)
  distce <- distce[genes.use, ]
  cell.set <- intersect(cell.set, colnames(distce))

  if (length(cell.set) > 1) {
    cv.vec <-  rowMeans(distce[names(genes.expr.prop), cell.set])
  } else {
    cv.vec <- distce[names(genes.expr.prop), cell.set]
  }
  names(cv.vec) <- names(genes.expr.prop)
  df <- data.frame(distance = cv.vec, expr.prop = genes.expr.prop)
  row.names(df) <- names(genes.expr.prop)
  return(df[order(df$distance), ])
}



## plot related functions------------------------------------

#' @importFrom ggplot2 ggsave scale_x_reverse scale_y_reverse
write_fig <- function(
  plt, filename = "myfigure", dir_name = "Figs", y_reverse = FALSE,
  x_reverse = FALSE, width = 7, height = 5.5, dpi = 200) {
  if (y_reverse) {
    plt <- plt + scale_y_reverse()
  }
  if (x_reverse) {
    plt <- plt + scale_x_reverse()
  }
  if (!is.null(filename)) {
    if (!dir.exists(dir_name)) {
      dir.create(dir_name)
    }
    ggplot2::ggsave(
      file = paste0("./", dir_name, "/", filename, ".png"), plot = plt,
      width = width, height = height, units = "in", dpi = dpi)
  }
}


#' Calculate UMAP projections for coembedding of cells and features
#' @description Calculate UMAP projections for coembedding of cells and features
#' @param seu a Seurat object with coembedding in the reductions slot wiht component name reduction.
#' @param reduction a string, specify the reduction component that denotes coembedding.
#' @param reduction.name a string, specify the reduction name for the obtained UMAP projection.
#' @param gene.set a string vector,  specify the features (genes) in calculating the UMAP projection, default as all features.
#' @param slot an optional string,  specify the slot in the assay, default as `data`.
#' @param assay an optional string,  specify the assay name in the Seurat object when adding the UMAP projection.
#' @param seed an optional integer, specify the random seed for reproducibility.
#' @return return a revised Seurat object by adding a new reduction component named `reduction.name`.
#' @details None
#' @seealso None
#' @references None
#' @export
#' @importFrom  Seurat Loadings Embeddings CreateDimReducObject
#'
#' @examples
#' data(pbmc3k_subset)
#' data(top5_signatures)
#' \donttest{
#' pbmc3k_subset <- coembedding_umap(
#'   pbmc3k_subset, reduction = "ncfm", reduction.name = "UMAPsig",
#'   gene.set = top5_signatures$gene
#' )
#' }
#'
#'
coembedding_umap <- function(seu, reduction, reduction.name, gene.set = NULL,
                             slot = "data", assay = "RNA", seed = 1) {
  
  
  calculateUMAP <- function(...){
    if (requireNamespace("scater", quietly = TRUE)) {
      x <- scater::calculateUMAP(...)
      return(x)
    } else {
      stop("coembedding_umap: scater is not available. Install scater to use this functionality.")
    } 
  }
  
  if (is.null(gene.set)) gene.set <- rownames(seu)
  gene.set <- unique(gene.set)
  if(!inherits(gene.set, "character"))
    stop("coembedding_umap: gene.set must be a character stering!")
  febd <- Loadings(seu, reduction)[gene.set, ]
  cebd <- Embeddings(seu, reduction)
  set.seed(seed)
  umap_all <- calculateUMAP(t(rbind(febd, cebd)))
  colnames(umap_all) <- paste0(gsub("_", "", reduction.name), "_", 1:2)
  n_gene <- length(gene.set)

  seu@reductions[[reduction.name]] <- CreateDimReducObject(
    embeddings = umap_all[-c(1:n_gene), ],
    loadings = umap_all[1:n_gene, ],
    key = paste0(gsub("_", "", reduction.name), "_"), assay = assay)
  seu
}

#' Coembedding dimensional reduction plot
#' @description Graph  output of a dimensional reduction technique on a 2D scatter plot where each point is a cell or feature and it's positioned based on the  coembeddings determined by the reduction technique. By default, cells and their signature features are colored by their identity class (can be changed with the group.by parameter).
#' @param seu a Seurat object with coembedding in the reductions slot wiht component name reduction.
#' @param reduction a string, specify the reduction component that denotes coembedding.
#' @param gene_txtdata a data.frame object with columns indcluding `gene` and `label`, specify the cell type/spatial domain and signature genes. Default as NULL, all features will be used in comebeddings.
#' @param cell_label an optional character in columns of metadata,  specify the group of cells/spots. Default as NULL, use Idents as the group.
#' @param xy_name an optional character,  specify the names of x and y-axis, default as the same as reduction.
#' @param dims a postive integer vector with length 2, specify the two components for visualization.
#' @param cols an optional string vector, specify the colors for cell group in visualization.
#' @param shape_cg a positive integers with length 2, specify the shapes of cell/spot and feature in plot.
#' @param pt_size an optional integer, specify the point size, default as 1.
#' @param pt_text_size an optional integer, specify the point size of text, default as 5.
#' @param base_size an optional integer, specify the basic size.
#' @param base_family an optional character, specify the font.
#' @param legend.point.size an optional integer, specify the point size of legend.
#' @param legend.key.size an optional integer, specify the size of legend key.
#' @param alpha an optional positive real, range from 0 to 1, specify the transparancy of points.
#' @return return a ggplot object
#' @details None
#' @seealso \code{\link{coembedding_umap}}
#' @references None
#' @export
#' @importFrom  ggplot2 ggplot aes_string geom_point scale_colour_manual scale_shape_manual theme_classic theme guides
#' @importFrom  Seurat Loadings  Embeddings Idents
#' @examples
#' data(pbmc3k_subset)
#' data(top5_signatures)
#' coembed_plot(pbmc3k_subset, reduction = "UMAPsig",
#'  gene_txtdata = top5_signatures,  pt_text_size = 3, alpha=0.3)
#'
coembed_plot <- function (seu, reduction, gene_txtdata = NULL, cell_label = NULL,
                          xy_name = reduction, dims = c(1, 2), cols = NULL, shape_cg=c(1,5),
                          pt_size=1,pt_text_size=5, base_size=16, base_family='serif',
                          legend.point.size=5, legend.key.size=1.5,alpha=0.3){
  # require(dplyr)
  # require(ggplot2)
  # require(ggrepel)
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    if (requireNamespace("grDevices", quietly = TRUE)) {
      x <- grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
      return(x)
    } else {
      stop("coembed_plot: grDevices is not available. Install grDevices to use plotting functionalities.")
    } 
  }
  geom_text_repel_here <- function(...){
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      x <- ggrepel::geom_text_repel(...)
      return(x)
    } else {
      stop("coembed_plot: ggrepel is not available. Install ggrepel to use plotting functionalities.")
    } 
  }
  if (is.null(gene_txtdata)) {
    umap_febd <- Loadings(seu, reduction)[, dims]
  }else {
    if (!all(c("gene", "label") %in% colnames(gene_txtdata))) {
      stop("'gene' and 'label' must be columns in 'gene_txtdata'!")
    }
    umap_febd <- Loadings(seu, reduction)[, dims]
    umap_febd <- umap_febd[gene_txtdata$gene, ]
  }
  umap_cebd <- Embeddings(seu, reduction)[, dims]
  umap_co_ebd <- as.data.frame(rbind(umap_febd, umap_cebd))
  colnames(umap_co_ebd) <- paste0(xy_name, dims)
  n_gene <- nrow(umap_febd)
  cgtype <- factor(c(rep("gene", n_gene), rep("cell", nrow(umap_cebd))))
  umap_co_ebd$cgtype <- cgtype

  if (is.null(cell_label)) {
    y <- as.character(Idents(seu))
  }else {
    y <- as.character(seu@meta.data[, cell_label])
  }
  if (is.null(cols)) {
    col_clusters <- gg_color_hue(length(unique(y)) +
                                             1)
    col_clusters[1] <- "#808080"
  }else {
    col_clusters <- cols
  }
  legend_ord <- c("gene", sort(unique(y)))
  names(col_clusters) <- legend_ord
  umap_co_ebd$cluster <- factor(c(rep("gene", n_gene),
                                  y), levels = legend_ord)

  p1 <- ggplot(data=umap_co_ebd, aes_string(x=colnames(umap_co_ebd)[1],
                                            y=colnames(umap_co_ebd)[2],
                                            shape='cgtype',color='cluster'))+
    geom_point(data = subset(umap_co_ebd, cgtype != "gene"), size = pt_size, alpha=alpha) + # First surface
    geom_point(data = subset(umap_co_ebd, cgtype == "gene"), size = pt_size+1.5) +
    scale_colour_manual(values=col_clusters) +
    scale_shape_manual(values = c("cell" = 16, "gene" = 4))+
    theme_classic(base_size=base_size, base_family = base_family)
  if (is.null(gene_txtdata)) {
    return(p1)
  }

  df11 <- data.frame(UMAP1 = umap_febd[, 1], UMAP2 = umap_febd[, 2],
                     label = row.names(umap_febd),
                     cgtype = rep("gene", nrow(umap_febd)),
                     cluster = rep("gene", nrow(umap_febd)), color = col_clusters[gene_txtdata$label])
  colnames(df11)[1:2] <-  paste0('xy_name', dims)
  p12 <- p1 + geom_text_repel_here(data = df11, aes_string(x=colnames(df11)[1],
                                                                 y=colnames(df11)[2],
                                                                 label = 'label'), size = pt_text_size,
                                         color = df11$color, max.overlaps = 50) +
      theme(legend.key.size = unit(legend.key.size, "lines"))+
      guides(shape = guide_legend(override.aes = list(size = legend.point.size)),
             color= guide_legend(override.aes = list(size = legend.point.size)))
  
  
  return(p12)
}


