#' @importFrom mclust Mclust mclustBIC 
mycluster <- function(Z, G, int.model='EEE', verbose=FALSE){
  ## require(mclust)
  mclus2 <- mclust::Mclust(Z, G=G,modelNames =int.model ,verbose=verbose)
  return(mclus2)
}

parafun_int <- function(k, Z, Sigma_equal, Sigma_diag, seed=1, init.start=5, int.model='EEE', verbose=FALSE){
  
  
  loglik0 = -1e20
  for (i in 1:init.start) {
    set.seed(seed+(i-1)*10)
    mclus0 <- mycluster(Z, G=k, int.model =int.model ,verbose=verbose)
    if (mclus0$loglik > loglik0) {
      loglik0 <- mclus0$loglik
      mclus2 <- mclus0
    }
  }
  
  Mu0k <- t(mclus2$parameters$mean)
  Sigmak <- mclus2$parameters$variance$sigma
  if(Sigma_diag){
    Sigma0k <- array(0, dim=dim(Sigmak))
    for(kk in 1:k){
      diag(Sigma0k[,,kk]) <- diag(Sigmak[,,kk])
    }
  } else {
    Sigma0k <- Sigmak
  }
  if(Sigma_equal){
    for(kk in 1:k){
      Sigma0k[,,kk] <- apply(Sigma0k, c(1,2), mean)
    }
  }
  Pi0k <- mclus2$parameters$pro
  return(list(y0k = mclus2$classification, Mu0k=Mu0k, Sigma0k=Sigma0k, Pi0k=Pi0k))
}

#' @importFrom parallel detectCores makeCluster clusterExport parLapply stopCluster
get_initial_value4seqK <- function(
  Kseq, hZ, Sigma_equal, Sigma_diag, seed=1, init.start=5, int.model='EEE', verbose=FALSE, coreNum){
  nK = length(Kseq)
  ## require(mclust)
  if (nK>1 & coreNum>1) {
    ## at most use 80% CPU
    cores <- min(c(nK, 10, parallel::detectCores()*0.8, coreNum))
    if (Sys.info()[1]=="Windows") {
      cl <- parallel::makeCluster(cores) # memory can not be shared type in Windows.
    } else {
      cl <- parallel::makeCluster(cores, type='FORK') # memory-shared type in linux or Mac.
    }
    
    message("Starting parallel computing initial values...")
    parallel::clusterExport(cl, varlist = c("Mclust", "mclustBIC"), envir=environment())
    # Run
    intList <- parallel::parLapply(
      cl, X=Kseq, parafun_int, Z=hZ, Sigma_equal=Sigma_equal, Sigma_diag=Sigma_diag, seed=seed,
      init.start=init.start, int.model=int.model, verbose=verbose)
    parallel::stopCluster(cl)
  } else {
    intList <- list(
      parafun_int(Kseq, Z=hZ, Sigma_equal=Sigma_equal, Sigma_diag=Sigma_diag, seed=seed,
                  init.start=init.start, int.model=int.model, verbose=verbose))
  }
  
  Mu0k     = lapply(intList, function(intListK) intListK$Mu0k)
  Sigma0k  = lapply(intList, function(intListK) intListK$Sigma0k)
  y0k      = lapply(intList, function(intListK) intListK$y0k)
  Pi0k     = lapply(intList, function(intListK) intListK$Pi0k)
  
  return(list(y0k = y0k, Mu0k=Mu0k, Sigma0k=Sigma0k, Pi0k=Pi0k))
}

#' @importFrom methods setClass new as
setClass("iSCMEBResObj", slots=list(
  posList = "ANY", 
  paramList= "list", 
  fitList = "ANY",
  project = "character",
  reduction = "list",
  idents = "ANY"
) )



fit.iscmeb <- function(
  VList, AdjList, K, beta_grid=seq(0, 5, by=0.2), maxIter_ICM=6, maxIter=25, 
  epsLogLik=1e-5, verbose=TRUE, int.model="EEE", init.start=1, 
  Sigma_equal = FALSE, Sigma_diag = TRUE, seed=1, coreNum=1, 
  criteria=c("MBIC", "BIC", "AIC"), c_penalty=1, dr.method="PCA"){
  
  degree.freedom <- function(K, q, nT, Sigma_equal, Sigma_diag, Sp_embed) {
    if (Sigma_diag) {
      # message("Sigma is set to diagonal matrices.\n")
      df_Sigma = q
    } else {
      # message("Sigma is set to dense matrices.\n")
      df_Sigma = q*(q+1)/2.0
    }
    if (Sigma_equal) {
      df_Sigma = df_Sigma
    } else {
      df_Sigma = df_Sigma*K
    }
    
    if (Sp_embed) {
      df_psi = q*(q+1)/2.0
    } else {
      df_psi = 0
    }
    df_psi = df_psi*nT
    
    # Mu + Sigma + Psi + beta
    dfree <- K*q + df_Sigma + df_psi + nT
    
    return(dfree)
  }
  
  error_heter <- TRUE
  M <- length(VList)
  q <- ncol(VList[[1]])
  K.order = order(K, decreasing = T)
  K = sort(K, decreasing = T)
  
  Psi_int <- array(dim=c(q,q, M))
  for( j in 1: M) Psi_int[,,j] <- cov(VList[[j]]/4)
  
  message("Evaluate initial values...")
  intlist <- get_initial_value4seqK(
    K, matlist2mat(VList), Sigma_equal=Sigma_equal, Sigma_diag=Sigma_diag, seed=seed, 
    init.start=init.start, int.model=int.model, verbose=verbose, coreNum=coreNum) 
  
  Mu_int    <- intlist$Mu0k 
  Sigma_int <- intlist$Sigma0k
  y_int     <- intlist$y0k
  
  message("Fit SC-MEB2...")
  Sp_embed <- TRUE
  result <- iSCMEBCpp(
    VList, AdjList, y_int, Mu_int, Sigma_int, Psi_int, 1.5, 
    beta_grid, maxIter_ICM, maxIter, epsLogLik, verbose, 
    !error_heter, Sigma_equal, Sigma_diag, Sp_embed, max(K), min(K), coreNum) 
  
  result <- result[order(K.order)]
  K  <- K[order(K.order)]
  
  output <- new(
    Class = "iSCMEBResObj",
    posList = NULL,
    paramList= list(),
    fitList = result,
    project = "iSC.MEB",
    reduction = list(),
    idents = NULL
  )
  
  n  <- sum(sapply(VList, nrow))
  nT <- length(VList)
  
  output@paramList$dfs = sapply(K, function(K0) degree.freedom(K0, q, nT, Sigma_equal, Sigma_diag, Sp_embed))
  output@paramList$log_likes = unlist(lapply(result, function(fit) fit$loglik))
  #output@paramList$MAIC <- -2.0*output@paramList$log_likes + output@paramList$dfs*2*log(log(q+n))*c_penalty
  output@paramList$AIC  <- -2.0*output@paramList$log_likes + output@paramList$dfs*2
  output@paramList$MBIC <- -2.0*output@paramList$log_likes + output@paramList$dfs*log(n)*log(log(q+n))*c_penalty
  output@paramList$BIC  <- -2.0*output@paramList$log_likes + output@paramList$dfs*log(n)
  
  criteria = match.arg(criteria)
  Mycriteria <- switch(
    criteria, 
    #MAIC = output@paramList$MAIC,
    AIC  = output@paramList$AIC,
    MBIC = output@paramList$MBIC, 
    BIC  = output@paramList$BIC 
  )
  
  output@paramList$opt  = which.min(Mycriteria)
  output@paramList$optK = K[output@paramList$opt]
  if (is.null(names(VList))) {
    output@paramList$sample_name = paste0("Sample", c(1:nT))
  } else {
    output@paramList$sample_name = names(VList)
  }
  output@paramList$K  = K
  output@paramList$n  = sapply(VList, nrow)
  output@paramList$q  = q
  output@paramList$Sigma_diag  = Sigma_diag
  output@paramList$Sigma_equal = Sigma_equal
  output@paramList$Sp_embed    = Sp_embed
  output@paramList$nT = nT
  output@paramList$modelselect = switch(
    criteria, 
    #MAIC = paste0("MAIC_", c_penalty),
    AIC  = "AIC",
    MBIC = paste0("MBIC_", c_penalty), 
    BIC  = "BIC" 
  )
  output@paramList$dr.method <- dr.method
  
  output@reduction[[dr.method]] <- VList
  output@reduction$iSCMEB = output@fitList[[output@paramList$opt]]$hZ
  output@idents = lapply(output@fitList[[output@paramList$opt]]$cluster, as.vector)
  
  return(output)
}

#' Fit an iSC-MEB model using specified multi-section embeddings
#' @description  Integrate multiple SRT data based on the PRECASTObj by ProFAST and iSC-MEB model fitting.
#' @param VList a M-length list of embeddings. The i-th element is a ni * q matrtix, where ni is the number of spots of sample i, and q is the number of embeddings. We provide this interface for those users who would like to define the embeddings by themselves.
#' @param AdjList an M-length list of sparse matrices with class \code{dgCMatrix}, specify the adjacency matrix used for intrisic CAR model in ProFAST. We provide this interface for those users who would like to define the adjacency matrix by themselves.
#' @param K an integer, specify the number of clusters.
#' @param beta_grid an optional vector of positive value, the candidate set of the smoothing parameter to be searched by the grid-search optimization approach, defualt as a sequence starts from 0, ends with 5, increase by 0.2.
#' @param maxIter the maximum iteration of ICM-EM algorithm. The default is 25.
#' @param epsLogLik a string, the species, one of 'Human' and 'Mouse'.
#' @param verbose an optional intger, spcify the number of housekeeping genes to be selected.
#' @param int.model an optional string, specify which Gaussian mixture model is used in evaluting the initial values for iSC.MEB, default as "EEE"; and see \code{\link{Mclust}} for more models' names.
#' @param init.start an optional number of times to calculate the initial value (1 by default). When init.start is larger than 1, initial value will be determined by log likelihood of mclust results.
#' @param Sigma_equal an optional logical value, specify whether Sigmaks are equal, default as \code{FALSE}.
#' @param Sigma_diag an optional logical value, specify whether Sigmaks are diagonal matrices, default as \code{TRUE}.
#' @param seed an optional integer, the random seed in fitting iSC-MEB model.
#' @return returns a iSCMEBResObj object which contains all model results.
#' @export
#'
iscmeb_run <- function(VList, AdjList, K, beta_grid=seq(0, 5, by=0.2),  maxIter=25, 
                       epsLogLik=1e-5, verbose=TRUE, int.model="EEE", init.start=1, 
                       Sigma_equal = FALSE, Sigma_diag = TRUE, seed=1){
  
  if(!is.list(VList)) stop("iscmeb_run: VList must be a list!")
  if(!is.list(AdjList)) stop("iscmeb_run: AdjList must be a list!")
  if(length(K)>1 || K<1) stop("iscmeb_run: K must be a postive integer!")
  
  
  reslist_iscmeb <-  fit.iscmeb( VList, AdjList, K, beta_grid=beta_grid, maxIter_ICM=6, maxIter=maxIter, 
                                 epsLogLik=epsLogLik, verbose=verbose, int.model= int.model, init.start=init.start, 
                                 Sigma_equal = Sigma_equal, Sigma_diag = Sigma_diag, seed=seed, coreNum=1, 
                                 criteria='MBIC', c_penalty=1, dr.method="PCA")
  return(reslist_iscmeb)
}





# iSC-MEB for PRECAST Object ---------------------------------------------
get_sampleID <- function(XList){
  sampleID <- list()
  r_max <- length(XList)
  for(r in 1:r_max){
    sampleID[[r]] <- rep(r, nrow(XList[[r]]))
  }
  sampleID <- unlist(sampleID)
  return(sampleID)
}
## use Harmony and Louvain clustering to select the number of clusters
#' @importFrom Seurat CreateSeuratObject CreateDimReducObject FindNeighbors FindClusters
drLouvain <- function(hZ, resolution=0.8){
  ### Louvain cluster based on estimated integrative low-dimensional embeddings. 
  
  ## require(Seurat)
  n <- nrow(hZ); q <- ncol(hZ)
  row.names(hZ) <- paste0("spot", 1:n)
  colnames(hZ) <- paste0("gene", 1:q)
  seu <- CreateSeuratObject(counts= t(hZ), assay='RNA')
  DefaultAssay(seu) <- "RNA"
  pca1 <- CreateDimReducObject(embeddings = hZ, key = "PC")
  seu@reductions$"pca" <- pca1
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:q, assay = "RNA")
  seu <- FindClusters(seu, resolution = resolution)
  return(seu$seurat_clusters)
}

#' Embedding alignment and clustering based on the embeddings from ProFAST
#' @description  Embedding alignment and clustering using the Harmony and Louvain based on the ebmeddings from ProFAST as well as determining the number of clusters.
#' @param PRECASTObj a PRECASTObj object created by \code{\link{CreatePRECASTObject}}.
#' @param resolution 	an optional real, the value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @return Return a revised \code{PRECASTObj} object with slot \code{PRECASTObj@resList} added by a \code{Harmony} compoonent (including the aligned embeddings and embeddings of batch effects) and a \code{Louvain} component (including the clusters).
#' @export
#' @importFrom harmony HarmonyMatrix
#'

RunHarmonyLouvain <- function(PRECASTObj, resolution=0.5){
  # require(Seurat)
  #require(harmony)
  
  if(is.null(PRECASTObj@resList$ProFAST)) stop("RunLouvain: please run ProFAST before using SelectclustNumber!")
  
  seed <- PRECASTObj@parameterList$seed
  verbose <- PRECASTObj@parameterList$verbose
  if(verbose)
    message( "******","Use Harmony to remove batch in the embeddings from ProFAST...")
  tstart <- Sys.time()
  
  sampleID <- get_sampleID(PRECASTObj@resList$ProFAST$hV)
  nvec <- sapply(PRECASTObj@resList$ProFAST$hV, nrow)
  
  set.seed(seed)
  tic <- proc.time()
  hZ_harmony_profastP <- HarmonyMatrix(matlist2mat(PRECASTObj@resList$ProFAST$hV), meta_data = data.frame(sample = sampleID),
                                       vars_use = "sample", do_pca = F)
  toc <- proc.time()
  .logDiffTime(sprintf(paste0("%s Finish Harmony correction"), "*****"), t1 = tstart, verbose = verbose)
  
  if(verbose)
    message( "******","Use Louvain to cluster and determine the number of clusters...")
  tstart <- Sys.time()
  res_louvain_harmony_profastP <- drLouvain(hZ_harmony_profastP, resolution = resolution)
  ## Get the batch embed in Harmony
  lm1 <- lm(matlist2mat(PRECASTObj@resList$ProFAST$hV)~hZ_harmony_profastP)
  batchembed <- mat2list(residuals(lm1), nvec = nvec)
  rm(lm1)
  PRECASTObj@resList$Harmony <- list()
  PRECASTObj@resList$Harmony$harmonyembed <- mat2list(hZ_harmony_profastP, nvec=nvec)
  PRECASTObj@resList$Harmony$batchEmbed <- batchembed
  PRECASTObj@resList$Louvain <- list()
  PRECASTObj@resList$Louvain$cluster <- vec2list(res_louvain_harmony_profastP, nvec=nvec)
  .logDiffTime(sprintf(paste0("%s Finsh Louvain clustering and find the optimal number of clusters"), "*****"), t1 = tstart, verbose = verbose)
  
  
  PRECASTObj@parameterList$K_opt <- length(unique(res_louvain_harmony_profastP))
  return(PRECASTObj)
  
}



#' Fit an iSC-MEB model using the embeddings from ProFAST
#' @description  Fit an iSC-MEB model using the embeddings from ProFAST and the number of clusters obtained by Louvain.
#' @param PRECASTObj a PRECASTObj object created by \code{\link{CreatePRECASTObject}}.
#' @param ... other arguments passed to \code{\link{iscmeb_run}}.
#' @return Return a revised PRECASTObj object with an added component \code{iSCMEB} in the slot \code{PRECASTObj@resList} (including the aligned embeddings, clusters and posterior probability matrix of clusters).
#' @export
#'
RuniSCMEB <- function(PRECASTObj, ...){
  
  verbose <- PRECASTObj@parameterList$verbose
  if(verbose)
    message( "******","Perform embedding alignment and spatial clustering using iSC-MEB based on  the embeddings obtained by ProFAST...")
  tstart <- Sys.time()
  
  
  tic <- proc.time()
  reslist_iscmeb <- iscmeb_run(
    PRECASTObj@resList$ProFAST$hV,
    PRECASTObj@AdjList,
    K=PRECASTObj@parameterList$K_opt, ...)
  toc <- proc.time()
  PRECASTObj@resList$iSCMEB <- list()
  PRECASTObj@resList$iSCMEB$cluster <- reslist_iscmeb@idents
  PRECASTObj@resList$iSCMEB$alignedEmbed <- reslist_iscmeb@reduction$iSCMEB
  PRECASTObj@resList$iSCMEB$batchEmbed <- reslist_iscmeb@fitList[[1]]$hV 
  PRECASTObj@resList$iSCMEB$RList <- reslist_iscmeb@fitList[[1]]$Rf
  .logDiffTime(sprintf(paste0("%s Finish iSC-MEB fitting"), "*****"), t1 = tstart, verbose = verbose)
  
  
  return(PRECASTObj)
}

