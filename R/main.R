
# generate man files
# devtools::document()
# R CMD check --as-cran ProFAST_1.1.tar.gz
## usethis::use_data(dat_r2_mac)
# pkgdown::build_site()
# pkgdown::build_home()
# pkgdown::build_article("ProFASTsimu")
# pkgdown::build_article("ProFASTdlpfc2")
# Basic functions ---------------------------------------------------------
.logDiffTime <- function(main = "", t1 = NULL, verbose = TRUE, addHeader = FALSE,
                         t2 = Sys.time(), units = "mins", header = "*****",
                         tail = "elapsed.", precision = 3)
{
  
  # main = ""; t1 = NULL; verbose = TRUE; addHeader = FALSE;
  # t2 = Sys.time(); units = "mins"; header = "###########";
  # tail = "elapsed."; precision = 3
  if (verbose) {
    timeStamp <- tryCatch({
      dt <- abs(round(difftime(t2, t1, units = units),
                      precision))
      if (addHeader) {
        msg <- sprintf("%s\n%s : %s, %s %s %s\n%s",
                       header, Sys.time(), main, dt, units, tail,
                       header)
      }
      else {
        msg <- sprintf("%s : %s, %s %s %s", Sys.time(),
                       main, dt, units, tail)
      }
      if (verbose)
        message(msg)
    }, error = function(x) {
      if (verbose)
        message("Time Error : ", x)
    })
  }
  
  return(invisible(0))
}


.logTime <- function(main='', prefix='*****', versoe=TRUE){
  
  if(versoe){
    message(paste0(Sys.time()," : ", prefix," ",  main))
  }
  
  
}


Diag<- function(vec){
  q <- length(vec)
  if(q > 1){
    y <- diag(vec)
  }else{
    y <- matrix(vec, 1,1)
  }
  return(y)
}

#' @importFrom irlba irlba
approxPCA <- function(X, q){ ## speed the computation for initial values.
  
  #rrequire(irlba) 
  n <- nrow(X)
  svdX  <- irlba(A =X, nv = q)
  PCs <- svdX$u %*% Diag(svdX$d[1:q])
  loadings <- svdX$v
  dX <- PCs %*% t(loadings) - X
  Lam_vec <- colSums(dX^2)/n
  return(list(PCs = PCs, loadings = loadings, Lam_vec = Lam_vec))
}


matlist2mat <- function (XList) 
{
  r_max <- length(XList)
  X0 <- XList[[1]]
  if (r_max > 1) {
    for (r in 2:r_max) {
      X0 <- rbind(X0, XList[[r]])
    }
  }
  return(X0)
}
mat2list <- function(z_int, nvec){
  
  zList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){
    
    zList_int[[i]] <- z_int[istart: sum(nvec[1:i]), ]
    istart <- istart + nvec[i]
  }
  return(zList_int)
}
vec2list <- function(y_int, nvec){
  if(length(y_int) != sum(nvec)) stop("vec2list: Check the argument: nvec!")
  
  yList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){
    
    yList_int[[i]] <- y_int[istart: sum(nvec[1:i])]
    istart <- istart + nvec[i]
  }
  return(yList_int)
}

get_indexList <- function(alist){
  nsample <- length(alist)
  nr <- 0
  indexList <- list()
  for(i in 1:nsample){
    indexList[[i]] <- (nr+1):(nrow(alist[[i]] )+nr)
    nr <- nr + nrow(alist[[i]] )
  }
  return(indexList)
}
## Define functions
#' (Varitional) ICM-EM algorithm for implementing ProFAST model
#' @description (Varitional) ICM-EM algorithm for implementing ProFAST model
#' @param XList an M-length list consisting of multiple matrices with class \code{dgCMatrix} or \code{matrix} that specifies the count/log-count gene expression matrix for each data batch used for ProFAST model.
#' @param AdjList an M-length list of sparse matrices with class \code{dgCMatrix}, specify the adjacency matrix used for intrisic CAR model in ProFAST. We provide this interface for those users who would like to define the adjacency matrix by themselves.
#' @param q an optional integer, specify the number of low-dimensional embeddings to extract in ProFAST. Larger q means more information extracted.
#' @param fit.model an optional string, specify the version of ProFAST to be fitted. The Gaussian version models the log-count matrices while the Poisson verions models the count matrices; default as \code{gaussian} due to fastter computation.
#' @param AList an optional list with each component being a vector whose length is equal to the rows of component in \code{XList}, specify the normalization factor in ProFAST. The default is \code{NULL} that means the normalization factor equal to 1. 
#' @param maxIter the maximum iteration of ICM-EM algorithm. The default is 30.
#' @param epsLogLik an optional positive vlaue, tolerance of relative variation rate of the observed pseudo loglikelihood value, defualt as '1e-5'.
#' @param error_heter a logical value, whether use the heterogenous error for ProFAST model, default as \code{TRUE}. If \code{error.heter=FALSE}, then the homogenuous error is used.
#' @param Psi_diag a logical value, whether set the conditional covariance matrix of the intrisic CAR to diagonal, default as \code{FALSE}.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed a postive integer, the random seed to be set in initialization.
#' @param Vint_zero an optional logical value, specify whether the intial value of intrisic CAR component is set to zero; default as \code{FALSE}.
#' @return return a list including the following components: (1) hV: an M-length list consisting of spatial embeddings in ProFAST; (2) nu: the estimated intercept vector; (3) Psi: the estimated covariance matrix; (4) W: the estimated shared loading matrix; (5) Lam: the estimated covariance matrix of error term; (6): ELBO: the ELBO value when algorithm convergence; (7) ELBO_seq: the ELBO values for all itrations.
#' @details None
#' @seealso \code{\link{ProFAST_structure}}, \code{\link{ProFAST}}, \code{\link{model_set_ProFAST}}
#' @references None
#' @export
#' @useDynLib ProFAST, .registration = TRUE
#' @importFrom  Matrix sparseMatrix
#'
#'
ProFAST_run <- function(XList, AdjList, q = 15,  fit.model = c("gaussian", "poisson"),
                       AList=NULL, maxIter = 25,
                       epsLogLik = 1e-5,verbose=TRUE, seed=1,
                       error_heter=TRUE, Psi_diag=FALSE, Vint_zero=FALSE){
  # XList is a sparse matrix
  # q = 15; maxIter = 25; seed=1;
  # epsLogLik = 1e-5;verbose=TRUE; error_heter=TRUE; Sigma_diag=TRUE
  #require(Matrix)
  ### Check arguments
  if(!is.list(XList)) stop("ProFAST_run: XList must be a list!")
  if(!is.list(AdjList)) stop("ProFAST_run: AdjList must be a list!")
  if(!is.list(AdjList) || !is.null(AList)) stop("ProFAST_run: AList must be a list or NULL!")
  if(q<1) stop("ProFAST_run: q must be an integer greater than 0!")
  
  
  for(r in seq_along(XList)){
    if(is.matrix(XList[[r]])){
      tmpMat <- XList[[r]]
      XList[[r]] <- as(tmpMat, "sparseMatrix")
    }
  }
  M <- length(XList)
  fit.model <- match.arg(fit.model)
  if(fit.model=='poisson'){
    XList_log <- lapply(XList, function(x) log(1+x))
    if(is.null(AList)){
      AList <- list(); # Normalization factor in Poisson factor model
      for(r in 1:M){
        AList[[r]] <- rep(1, nrow(XList[[r]]))
      } 
    }
  }else if(fit.model == 'gaussian'){
    XList_log <- XList
  }
  
  XList_center <- lapply(XList_log, scale, scale=FALSE)
  rm(XList_log)
  nv_int <- t(sapply(XList_center, function(x) attr(x, "scaled:center")))
  
  set.seed(seed)
  princ <- approxPCA(matlist2mat(XList_center), q= q)
  
  Zmat <- princ$PCs
  
  indexList <- get_indexList(XList)
  Psi_int <- array(0, dim=c(q,q, length(XList)))
  for(r in 1:M){
    tmp_mat <- cov(Zmat[indexList[[r]],])
    if(Psi_diag){
      diag(Psi_int[,,r]) <- diag(tmp_mat)
    }else{
      Psi_int[,,r] <- tmp_mat
    }
    
  }
  W_int <- princ$loadings# loading list
  EvList <- mat2list(Zmat, nvec=sapply(XList, nrow))
  
  rm(princ)
  p <- ncol(XList[[1]])
  LamMat_int <- matrix(NA, nrow=M, ncol=p, byrow=T)
  for(r in 1: M){
    if(error_heter){
      LamMat_int[r, ] <- apply(XList_center[[r]] - EvList[[r]] %*% t(W_int), 2, var)
    }else{
      LamMat_int[r, ] <- apply(matlist2mat(XList_center) - Zmat %*% t(W_int), 2, var)
    }
  }
  rm(XList_center)
  if(Vint_zero){
    for(r in seq_along(EvList)){
      EvList[[r]] <- matrix(0, nrow(EvList[[r]]), q)
    }
  }
  
  if(fit.model=='poisson'){
    
    reslist <- profast_p_cpp(Xlist= XList, AList=AList, Adjlist=AdjList, nv_int=nv_int, W_int, 
                             Lam_int=LamMat_int, Psi_int, 
                             EvList=EvList, maxIter=maxIter, 
                             epsELBO=epsLogLik, verbose=verbose, homo = !error_heter, Psi_diag=Psi_diag)
    
  }else if(fit.model == 'gaussian'){
    reslist <- profast_g_cpp(Xlist=XList, Adjlist=AdjList, nv_int, W_int, Lam_int=LamMat_int, Psi_int, 
                             EvList=EvList, maxIter=maxIter, 
                             epsLogLik=epsLogLik, verbose, homo = !error_heter, Psi_diag=Psi_diag)
    
  }
  ## Put the intercept term to nu
  reslist$nu <- reslist$nu + t(reslist$W %*% sapply(reslist$hV, colMeans))
  reslist$hV <- lapply(reslist$hV, scale, scale=FALSE)
  
  return(reslist)
} 




# Metrics -----------------------------------------------------------------



#'  Calcuate the the adjusted McFadden's pseudo R-square 
#' @description  Calcuate the the adjusted McFadden's pseudo R-square  between the  embeddings and the  labels
#' @param embeds a n-by-q matrix, specify the embedding matrix.
#' @param y a n-length vector, specify the labels.
#' @return return the adjusted McFadden's pseudo R-square. 
#' @details None
#' @references McFadden, D. (1987). Regression-based specification tests for the multinomial logit model. Journal of econometrics, 34(1-2), 63-82.
#' @export
#' @importFrom  nnet multinom
#' @importFrom performance r2_mcfadden
#'
get_r2_mcfadden <- function(embeds, y){
  # library(nnet)
  # library(performance)
  dat_r2_mac <- NULL
  
  # data('dat_r2_mac', package = "ProFAST")
  y <- as.numeric(as.factor(y))
  model1 <- nnet::multinom(y~embeds)
  R2 <- r2_mcfadden(model1)
  return(R2$R2_adjusted)
}


# Design high-level function using PRECAST objects -------------------------
#' Set parameters for ProFAST model
#' @description  Prepare parameters setup for ProFAST model fitting.
#' @param maxIter the maximum iteration of ICM-EM algorithm. The default is 30.
#' @param epsLogLik an optional positive vlaue, tolerance of relative variation rate of the observed pseudo loglikelihood value, defualt as '1e-5'.
#' @param error_heter a logical value, whether use the heterogenous error for ProFAST model, default as \code{TRUE}. If \code{error.heter=FALSE}, then the homogenuous error is used.
#' @param Psi_diag a logical value, whether set the conditional covariance matrices of intrisic CAR to diagonal, default as \code{FALSE}
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed a postive integer, the random seed to be set in initialization.
#' @export
#' @examples
#' model_set_ProFAST(maxIter = 30, epsLogLik = 1e-5,
#'   error_heter=TRUE, Psi_diag=FALSE, verbose=TRUE, seed=2023)
#'
model_set_ProFAST <- function(maxIter = 30, epsLogLik = 1e-5,
                              error_heter=TRUE, Psi_diag=FALSE, verbose=TRUE, seed=1){
  
  para_settings <- list(maxIter = maxIter, seed=seed,
                        epsLogLik = epsLogLik,verbose=verbose,
                        error_heter=error_heter, Psi_diag=Psi_diag)
  return(para_settings)
}


#' (Varitional) ICM-EM algorithm for implementing ProFAST model with structurized parameters
#'
#' @param XList an M-length list consisting of multiple matrices with class dgCMatrix or matrix that specify the count/log-count gene expression matrix for each data batch used for ProFAST model.
#' @param AdjList an M-length list of sparse matrices with class dgCMatrix, specify the adjacency matrix used for intrisic CAR model in ProFAST. We provide this interface for those users who would like to define the adjacency matrix by themselves.
#' @param q an optional integer, specify the number of low-dimensional embeddings to extract in ProFAST
#' @param fit.model an optional string, specify the version of ProFAST to be fitted. The Gaussian version models the log-count matrices while the Poisson verions models the count matrices; default as gaussian due to fastter computation.
#' @param parameterList an optional list, specify other parameters in ProFAST model; see \code{\link{model_set_ProFAST}} for other paramters. The default is \code{NULL} that means the default parameters produced by \code{\link{model_set_ProFAST}} is used.
#' @return return a list including the following components: (1) hV: an M-length list consisting of spatial embeddings in ProFAST; (2) nu: the estimated intercept vector; (3) Psi: the estimated covariance matrix; (4) W: the estimated shared loading matrix; (5) Lam: the estimated covariance matrix of error term; (6): ELBO: the ELBO value when algorithm convergence; (7) ELBO_seq: the ELBO values for all itrations.
#' @details None
#' @seealso \code{\link{ProFAST_run}}, \code{\link{ProFAST}}, \code{\link{model_set_ProFAST}}
#' @references None
#' @export
#'
#'
ProFAST_structure <- function(XList, AdjList, q= 15,  fit.model = c("poisson", "gaussian"), 
                               parameterList = NULL){
  
  if(is.null(parameterList)){
    parameterList <- model_set_ProFAST()
  }
  ### initialize para: 6 arguments.
  maxIter <- epsLogLik<- verbose<- NULL
  error_heter<- Psi_diag<-  seed<-  NULL
  para_names <- c("maxIter", "epsLogLik", "verbose", "error_heter", "Psi_diag", "seed")
  
  for(iname in para_names){
    assign(iname, parameterList[[iname]])
  }
  fit.model <- match.arg(fit.model)
  
  
  resList <- ProFAST_run(XList=XList, AdjList=AdjList, q = q,  fit.model = fit.model,
                         AList= NULL, maxIter = maxIter, seed=seed,
                         epsLogLik = epsLogLik,verbose=verbose,
                         error_heter=error_heter, Psi_diag=Psi_diag, Vint_zero=FALSE)
  return(resList)
}

#' Add ProFAST model settings for a PRECASTObj object
#'
#' @param PRECASTObj a PRECASTObj object created by \code{\link{CreatePRECASTObject}}.
#' @param ... other arguments to be passed to \code{\link{model_set_ProFAST}} function.
#' @references None
#' @return  Return a revised PRECASTObj object with slot \code{parameterList} changed.
#' @export
#'
#'
AddParSettingProFAST <- function(PRECASTObj, ...){
  PRECASTObj@parameterList <- model_set_ProFAST(...)
  return(PRECASTObj)
}
  
#' Run ProFAST model for a PRECASTObj object
#'
#' @param PRECASTObj a PRECASTObj object created by \code{\link{CreatePRECASTObject}}.
#' @param q an optional integer, specify the number of low-dimensional embeddings to extract in ProFAST
#' @param fit.model an optional string, specify the version of ProFAST to be fitted. The Gaussian version models the log-count matrices while the Poisson verions models the count matrices; default as poisson.
#' @references None
#' @return  Return a revised PRECASTObj object with slot \code{PRECASTObj@resList} added by a \code{ProFAST} compoonent.
#' @export
#' @importFrom Matrix t
#' @importFrom Seurat DefaultAssay 
#'
#'

ProFAST <- function(PRECASTObj, q= 15, fit.model=c("poisson", "gaussian")){
  # suppressMessages(rrequire(Matrix))
  # suppressMessages(rrequire(Seurat))
  
  if(!inherits(PRECASTObj, "PRECASTObj")) 
    stop("ProFAST: Check the argument: PRECASTObj!  PRECASTObj must be a PRECASTObj object.")
  
  if(q < 0) stop("ProFAST: Check the argument: q!  PRECASTObj must be a positive integer.")
  
  if(is.null(PRECASTObj@seulist)) stop("The slot seulist in PRECASTObj is NULL!")
  
  fit.model <- match.arg(fit.model)
  
  verbose <- PRECASTObj@parameterList$verbose
  if(verbose){
    if(fit.model=="poisson"){
      message( "******","Run the Poisson version of ProFAST...")
    }else if(fit.model== "gaussian"){
      message( "******","Run the Gaussian version of ProFAST...")
    }else{
      stop("ProFAST: Check the argument: fit.model! It is not supported for this fit.model!")
    }
  }
    
  tstart <- Sys.time()
  ## Get normalized data 
  get_data <- function(seu, assay = NULL, fit.model='poisson'){
    if(is.null(assay)) assay <- DefaultAssay(seu)
    if(fit.model=='poisson'){
      dat <- Matrix::t(seu[[assay]]@counts)
    }else if(fit.model=='gaussian'){
      dat <- Matrix::t(seu[[assay]]@data)
    }else{
      stop("ProFAST: Check the argument: fit.model! It is not supported for this fit.model!")
    }
    
    return(dat)
  }
  XList <- lapply(PRECASTObj@seulist,  get_data, fit.model=fit.model)
  
  # PRECASTObj@resList$ProFAST <- list()
  PRECASTObj@resList$ProFAST <- ProFAST_structure(XList, q= 15,  fit.model = fit.model, 
                                          AdjList = PRECASTObj@AdjList, parameterList = PRECASTObj@parameterList)
  
  .logDiffTime(sprintf(paste0("%s Finish ProFAST"), "*****"), t1 = tstart, verbose = verbose)
  return(PRECASTObj)
}





# Select housekeeping genes for preparation of removing unwanted--------
#' Transfer gene names from one fortmat to the other format
#' @description  Transfer gene names from one fortmat to the other format for two species: human and mouse.
#' @param genelist a string vector, the gene list to be transferred.
#' @param now_name a string, the current format of gene names, one of 'ensembl', 'symbol'.
#' @param to_name a string, the  format of gene names to transfer, one of 'ensembl', 'symbol'.
#' @param species a string, the species, one of 'Human' and 'Mouse'.
#' @param Method a string, the method to use, one of 'biomaRt' and 'eg.db', default as 'eg.db'.
#' @return Return a string vector of transferred gene names. The gene names not matched in the database will not change.
#' @export
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom AnnotationDbi mapIds
#' @importFrom biomaRt useDataset getBM useMart
#' @examples
#' geneNames <- c("ENSG00000171885", "ENSG00000115756")
#' transferGeneNames(geneNames, now_name = "ensembl", to_name="symbol",species="Human", Method='eg.db')
#'
#'
transferGeneNames <- function(genelist, now_name = "ensembl", to_name="symbol",
                              species=c("Human", "Mouse"), Method=c('eg.db', 'biomart') ){
  
  firstup <- function(x) {
    ## First letter use upper capital
    x <- tolower(x)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  species <- match.arg(species)
  Method <- match.arg(Method)
  if(! toupper(species) %in% c("HUMAN", "MOUSE")) stop("Check species: the current version only support Human and Mouse!")
  transferredNames <- switch (toupper(species),
                              HUMAN = {
                                if(tolower(Method)=='eg.db'){
                                  #require(org.Hs.eg.db)
                                  mapIds(org.Hs.eg.db, keys = genelist,
                                         keytype = toupper(now_name), column=toupper(to_name))
                                }else if(tolower(Method)=='biomart'){
                                  #require(biomaRt)
                                  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
                                  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                                                  values=genelist,mart= mart)
                                  
                                  idx_in_genelst <- which(G_list$ensembl_gene_id %in% genelist)
                                  G_list <- G_list[idx_in_genelst,]
                                  idx_dup <- which(duplicated(G_list$ensembl_gene_id))
                                  G_list <- G_list[-idx_dup,]
                                  row.names(G_list) <- G_list$ensembl_gene_id
                                  symbol_list <- G_list[genelist,]$hgnc_symbol
                                  
                                  symbol_list[symbol_list==''] <- NA
                                  symbol_list
                                  
                                }else{
                                  stop("Check Method: the current version only support biomaRt and eg.db!")
                                }
                                
                                
                              },
                              MOUSE= {
                                if(tolower(Method)=='eg.db'){
                                  #rrequire(org.Mm.eg.db)
                                  mapIds(org.Mm.eg.db, keys = genelist,
                                         keytype = toupper(now_name), column=toupper(to_name))
                                }else if(tolower(Method)=='biomart'){
                                  #rrequire(biomaRt)
                                  mart <- useDataset(" mmusculus_gene_ensembl", useMart("ensembl"))
                                  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                                                  values=genelist,mart= mart)
                                  
                                  idx_in_genelst <- which(G_list$ensembl_gene_id %in% genelist)
                                  G_list <- G_list[idx_in_genelst,]
                                  idx_dup <- which(duplicated(G_list$ensembl_gene_id))
                                  G_list <- G_list[-idx_dup,]
                                  row.names(G_list) <- G_list$ensembl_gene_id
                                  symbol_list <- G_list[genelist,]$hgnc_symbol
                                  
                                  symbol_list[symbol_list==''] <- NA
                                  symbol_list
                                }else{
                                  stop("Check Method: the current version only support biomaRt and eg.db!")
                                }
                                
                              }
  )
  
  if(toupper(to_name) == 'SYMBOL'){
    if(toupper(species) == "MOUSE"){
      transferredNames <- firstup(transferredNames)
    }else{
      transferredNames <- toupper(transferredNames)
    }
    
  }
    
  
  flag_na <- is.na(transferredNames)
  if(any(flag_na))
    transferredNames[flag_na] <- genelist[flag_na]
  
  
  return(transferredNames)
}




selectHKFeatures <- function(seulist, HKFeatureList, HKFeatures=200){
  ## This function is used for selecting common informative features
  if(length(seulist) != length(HKFeatureList)) stop("The length of suelist and HKFeatureList must be equal!")
  if(length(seulist) ==1){
    if(length(HKFeatureList[[1]]) >= HKFeatures){
      genelist <- HKFeatureList[[1]][1:HKFeatures]
    }else{
      genelist <- HKFeatureList[[1]]
      warning("The IntFeature is larger than the  number of elements in FeatureList!")
    }
    return(genelist)
  } 
  geneUnion <- unique(unlist(HKFeatureList))
  ## ensure each seuobject has the genes in geneUnion
  
  # Remove zero-variance genes
  genes_zeroVar <- unique(unlist(lapply(seulist, function(x){
    assay <- DefaultAssay(x)
    geneUnion[Matrix::rowSums(x[[assay]]@counts[geneUnion,])==0]
  })))
  
  
  #geneUnion[pbapply::pbapply(x@assays$RNA@counts[geneUnion,],1, sd)==0])))
  gene_Var <- setdiff(geneUnion, genes_zeroVar)
  
  # sort by number of datasets that identified this gene as the gene without spatial variation
  nsample <- length(seulist)
  numVec <- rep(0, length(gene_Var))
  rankMat <-matrix(NA,length(gene_Var), nsample)
  row.names(rankMat) <- gene_Var
  for(i in 1:length(gene_Var)){
    for(j in 1:nsample){
      if(is.element(gene_Var[i], HKFeatureList[[j]])){
        numVec[i] <- numVec[i] +1
        rank1 <- which(HKFeatureList[[j]]==gene_Var[i])
        rankMat[i, j] <- rank1
      }
    }
    
  }
  
  cutNum <- sort(numVec, decreasing = T)[min(HKFeatures, length(numVec))]
  if(max(numVec)> cutNum){
    genelist1 <- gene_Var[numVec>cutNum]
  }else{
    genelist1 <- NULL
  }
  num_rest_genes <- min(HKFeatures, length(numVec)) - length(genelist1)
  
  gene2 <- gene_Var[numVec==cutNum]
  
  
  rankMat2 <- rankMat[gene2, ]
  rowMedian <- function(xmat, na.rm=TRUE){
    apply(xmat, 1, median, na.rm=na.rm)
  }
  genes1rank <- gene2[order(rowMedian(rankMat2, na.rm=T))[1:num_rest_genes]]
  genelist <- c(genelist1, genes1rank)
  
  return(genelist)
}


# Select housekeeping genes for preparation of removing unwanted--------
#' Select housekeeping genes
#' @description  Select housekeeping genes for preparation of removing unwanted variations in expression matrices
#' @param seuList an M-length  list consisting of Seurat object, include the information of expression matrix and spatial coordinates (named \code{row} and \code{col}) in the slot \code{meta.data}.
#' @param species a string, the species, one of 'Human' and 'Mouse'.
#' @param HK.number an optional integer, specify the number of housekeeping genes to be selected.
#' @return Return a string vector of the selected gene names. 
#' @export
#' @importFrom  DR.SC FindSVGs
#' @importFrom utils data
#' @importFrom pbapply pblapply
#'
SelectHKgenes <- function(seuList, species= c("Human", "Mouse"), HK.number=200){
  
  #rrequire(PRECAST)
  
  
  gene_symbols <- Reduce(intersect, lapply(seuList, row.names))
  # gene_symbols <- transferGeneNames(gene_intersect, species="Human")
  
  if(tolower(species) =="human"){
    #data("Human_HK_genes", package = "PRECAST")
    idx_hk <- which(gene_symbols %in% as.character(PRECAST::Human_HK_genes$Gene))
    housekeep_genes <- gene_symbols[idx_hk]
    
  }else if(tolower(species) == "mouse"){
    # data("Mouse_HK_genes", package = "PRECAST")
    idx_hk <- which(gene_symbols %in% as.character(PRECAST::Mouse_HK_genes$Gene))
    housekeep_genes <- gene_symbols[idx_hk]
  }
  
  if(length(housekeep_genes) < HK.number) HK.number <- length(housekeep_genes)
  if(length(housekeep_genes) < 5) warning("SelectHKgenes: the housekeeping genes are less than 5!")
  if(length(housekeep_genes) < 100){
    return(housekeep_genes)
  }
  
  seulist_HK <- lapply(seuList, function(x) x[housekeep_genes,])
  rm(seuList)
  seulist_HK <- pbapply::pblapply(seulist_HK, DR.SC::FindSVGs, nfeatures = nrow(seulist_HK[[1]]))
  geneList_noSpa <- lapply(seulist_HK, function(x) {
    
    order_idx <- order(x[['RNA']]@meta.features$adjusted.pval.SVGs, decreasing = T) # rank the gene with largest p-value the first
    genes <- row.names(x[['RNA']]@meta.features)[order_idx[1:HK.number]]
    # idx <- which(x[['RNA']]@meta.features$adjusted.pval.SVGs[order_idx]>0.05)
    #return(genes[idx])
    return(genes)
  }) 
  
  final_HKgenes <- selectHKFeatures(seulist_HK, geneList_noSpa, HKFeatures=HK.number)
  return(final_HKgenes)
}


# Integration analysis ----------------------------------------------------
#' @importFrom stats as.formula  coef cov lm median model.matrix model.matrix.lm residuals var 
correct_genesR <- function(XList, RList, HList, Tm, AdjList, covariateList=NULL, maxIter=30, epsELBO=1e-4, verbose=TRUE){
  
  dfList2df <- function(dfList){
    
    df <- dfList[[1]]
    r_max <- length(dfList)
    if(r_max>1){
      for(r in 2:r_max){
        df <- rbind(df, dfList[[r]])
      }
    }
    
    return(df)
  }
 
  M <- length(XList);  p <- ncol(XList[[1]])
  q <- ncol(HList[[1]]); d <- ncol(Tm)
  Xmat <- matlist2mat(XList)
  R <- matlist2mat(RList)
  nvec <- sapply(XList, nrow)
  
  if(!is.null(covariateList)){
    # covariates <- matlist2mat(covariateList)
    # covariates <- as.matrix(covariates)
    covarites_df <- dfList2df(covariateList)
    covariates <-  model.matrix.lm(object = ~.+1, data = covarites_df, na.action = "na.pass")
    rm(covariateList, covarites_df)
    R <- cbind(R, covariates[,-1])
    
    RList <- mat2list(R, nvec = nvec)
    rm(covariates)
  }
  colnames(R) <- NULL
  
  K <- ncol(RList[[1]]); 
  H <- matlist2mat(HList)
  colnames(H) <- NULL
  nvec <- sapply(RList, nrow)
  TmList <- list()
  for(m in 1:M){
    TmList[[m]] <- matrix(Tm[m,], nrow=nvec[m], ncol=ncol(Tm), byrow=T)
  }
  TM <- matlist2mat(TmList)
  
  lm1 <- lm(Xmat~ R+H+TM+0)
  coef_all <- coef(lm1)
  rm(R, H, Xmat, lm1)
  alphaj_int <- coef_all[paste0("R", 1:K),]
  gammaj_int <- coef_all[paste0("H", 1:q),]
  if(d == 1){
    zetaj_int <- matrix(coef_all[paste0("TM"), ], nrow=1)
  }else{
    zetaj_int <- coef_all[paste0("TM", 1:d)]
  }
  if(sum(Tm)<1e-20){
    if(ncol(Tm)>1){
      stop("Tm must be a one-column matrix when it is full-zero!")
    }else{
      zetaj_int <- matrix(1, 1, p)
    }
  }
  
  
  sigmaj_int <- matrix(1, M, p)
  psij_int <- matrix(1,M, p)
  # maxIter <- 30; epsELBO <- 1e-4; verbose<- TRUE
  reslist <- correct_genes(XList, RList, HList, Tm, Adjlist=AdjList, sigmaj_int, psij_int, 
                           alphaj_int, gammaj_int, zetaj_int, maxIter, epsELBO, verbose) 
  # str(reslist)
  correct_list <- list(XList_correct=lapply(1:M, function(j) XList[[j]] - HList[[j]] %*% reslist$gamma),
                       gammaj=reslist$gamma, alphaj=reslist$alpha, zetaj=reslist$zeta)
  

  
  return(correct_list)
  
}
correct_genes_subsampleR <- function(XList, RList, HList, Tm, AdjList, subsample_rate, 
                                     covariateList=NULL, maxIter=30, epsELBO=1e-4, verbose=TRUE){
  
  dfList2df <- function(dfList){
    
    df <- dfList[[1]]
    r_max <- length(dfList)
    if(r_max>1){
      for(r in 2:r_max){
        df <- rbind(df, dfList[[r]])
      }
    }
    
    return(df)
  }
  
  M <- length(XList);  p <- ncol(XList[[1]])
  K <- ncol(RList[[1]]); q <- ncol(HList[[1]]); d <- ncol(Tm)
  Xmat <- matlist2mat(XList)
  R <- matlist2mat(RList)
  colnames(R) <- NULL
  
  if(!is.null(covariateList)){
    # covariates <- matlist2mat(covariateList)
    # covariates <- as.matrix(covariates)
    covarites_df <- dfList2df(covariateList)
    covariates <-  model.matrix.lm(object = ~.+1, data = covarites_df, na.action = "na.pass")
    rm(covariateList, covarites_df)
    R <- cbind(R, covariates[,-1])
    rm(covariates)
  }
  
  # subsample_rate <- 0.1
  index_List <- get_indexList(RList)
  set.seed(1)
  index_subsample <- sort(sample(sum(nvec), sum(nvec)*subsample_rate))
  ## calculate the number of indices belonging to the index of each slide
  nvec_subsample <- rep(NA, length(nvec))
  for(i in 1:length(nvec_subsample)){
    message("i = ", i)
    
    nvec_subsample[i] <- sum(index_subsample%in% index_List[[i]])
  }
  index_List_new <- lapply(RList, function(x) 1: nrow(x))
  index_subsample_new <- unlist(index_List_new)[index_subsample]
  index_subsampleList <- vec2list(index_subsample_new, nvec_subsample)
  
  XList_sub <- list(); RList_sub <- list(); HList_sub <- list()
  AdjList_sub <- list()
  for(i in 1:length(XList)){
    message("i = ", i)
    index_tmp <- index_subsampleList[[i]]
    XList_sub[[i]] <- XList[[i]][index_tmp,]
    RList_sub[[i]] <- RList[[i]][index_tmp, ]
    HList_sub[[i]] <- HList[[i]][index_tmp, ]
    AdjList_sub[[i]] <- AdjList[[i]][index_tmp, index_tmp]
  }
  
  H <- matlist2mat(HList)
  colnames(H) <- NULL
  nvec <- sapply(RList, nrow)
  TmList <- list()
  for(m in 1:M){
    TmList[[m]] <- matrix(Tm[m,], nrow=nvec[m], ncol=ncol(Tm), byrow=T)
  }
  TM <- matlist2mat(TmList)
  
  lm1 <- lm(Xmat~ R+H+TM+0)
  coef_all <- coef(lm1)
  rm(R, H, Xmat, lm1)
  alphaj_int <- coef_all[paste0("R", 1:K),]
  gammaj_int <- coef_all[paste0("H", 1:q),]
  if(d == 1){
    zetaj_int <- matrix(coef_all[paste0("TM"), ], nrow=1)
  }else{
    zetaj_int <- coef_all[paste0("TM", 1:d)]
  }
  if(sum(Tm)<1e-20){
    if(ncol(Tm)>1){
      stop("Tm must be a one-column matrix when it is full-zero!")
    }else{
      zetaj_int <- matrix(1, 1, p)
    }
  }
  
  
  sigmaj_int <- matrix(1, M, p)
  psij_int <- matrix(1,M, p)
  
  
  # maxIter <- 30; epsELBO <- 1e-4; verbose<- TRUE
  reslist <- correct_genes(XList, RList, HList, Tm, Adjlist=AdjList, sigmaj_int, psij_int, 
                           alphaj_int, gammaj_int, zetaj_int, maxIter, epsELBO, verbose) 
  
  # str(reslist)
  correct_list <- list(XList_correct=lapply(1:M, function(j) XList[[j]] - HList[[j]] %*% reslist$gamma),
                       gammaj=reslist$gamma, alphaj=reslist$alpha, zetaj=reslist$zeta)
  
  
  # for(m in 1:M){
  #   if(!any(sapply(XList, function(x) is.null(row.names(x))))){
  #     
  #     row.names( correct_list$XList_correct[[m]]) <- paste0("slice", m, "_", row.names(XList[[m]]) )
  #   }
  #   if(!any(sapply(XList, function(x) is.null(colnames(x))))){
  #        colnames(correct_list$XList_correct[[m]]) <- colnames(XList[[m]])
  #   }
  # }
  
  return(correct_list)
  
}



#' Integrate multiple SRT data into a Seurat object
#' @description  Integrate multiple SRT data based on the \code{PRECASTObj} object by ProFAST and other model fitting.
#' @param PRECASTObj a PRECASTObj object created by \code{\link{CreatePRECASTObject}}.
#' @param seulist_HK a list with Seurat object as component including only the housekeeping genes.
#' @param Method a string, specify the method to be used and two methods are supprted: \code{iSC-MEB} and \code{HarmonyLouvain}. The default is \code{iSC-MEB}.
#' @param seuList_raw an optional list with Seurat object, the raw data.
#' @param covariates_use a string vector, the colnames in \code{PRECASTObj@seulist[[1]]@meta.data}, representing other biological covariates to considered when removing batch effects. This is achieved by adding additional covariates for biological conditions in the regression, such as case or control. Default as 'NULL', denoting no other covariates to be considered.
#' @param Tm an optional numeric vector with the length equal to \code{PRECASTObj@seulist}, the time point information if the data include the temporal information. Default as \code{NULL} that means there is no temporal information.
#' @param verbose an optional logical value, default as \code{TRUE}.
#' @return Return a Seurat object by integrating all SRT data batches into a SRT data, where the column "batch" in the meta.data represents the batch ID, and the column "cluster" represents the clusters. The embeddings are put in \code{seu@reductions} slot and \code{Idents(seu)} is set to cluster label. Note that only the normalized expression is valid in the data slot while count is invalid.
#' @export
#' @details If \code{seuList_raw} is not equal \code{NULL} or \code{PRECASTObj@seuList} is not \code{NULL}, this function will remove the unwanted variations for all genes in \code{seuList_raw} object. Otherwise, only the the unwanted variation of genes in \code{PRECASTObj@seulist} will be removed. The former requies a big memory to be run, while the latter note.
#' @importFrom Matrix t sparseMatrix
#' @importFrom Seurat DefaultAssay CreateSeuratObject `DefaultAssay<-` `Idents<-`
#' @importFrom PRECAST Add_embed
#' @import gtools
#' @useDynLib ProFAST, .registration = TRUE
#'
IntegrateSRTData <- function(PRECASTObj, seulist_HK, Method=c("iSC-MEB", "HarmonyLouvain"), seuList_raw=NULL, covariates_use=NULL, 
                             Tm=NULL, verbose=TRUE){
  # require(Matrix)
  # library(Seurat)
  # require(PRECAST)

  Method <- match.arg(Method)
  verbose <- verbose
  if(verbose)
    message( "******","Perform PCA on housekeeping gene expression matrix...")
  tstart <- Sys.time()
  nvec <- sapply(PRECASTObj@seulist, ncol)
  if(!is.null(seulist_HK)){
    if(Method=='iSC-MEB'){
      message("Use housekeeping genes and the results in iSC-MEB to remove unwanted variations...")
    }else if(Method=='HarmonyLouvain'){
      message("Use housekeeping genes and the results in Harmony and Louvain to remove unwanted variations...")
    }else{
      stop("IntegrateSRTData: it does not support this Method!")
    }
    M <- length(seulist_HK)
    seulist_HK <- lapply(1:M, function(r){
      seu_tmp <- seulist_HK[[r]]
      seu_tmp <- seu_tmp[, colnames(PRECASTObj@seulist[[r]])]
      return(seu_tmp)
    })
    XList_hk <- pbapply::pblapply(seulist_HK, function(x) Matrix::t(x[['RNA']]@data))
    nvec <- sapply(XList_hk, nrow)
    XList_hk <- pbapply::pblapply(XList_hk, as.matrix)
    Xmat_hk <- matlist2mat(XList_hk)
    princ <- approxPCA(Xmat_hk, q=10)
    HList <- mat2list(princ$PCs, nvec)
    rm(Xmat_hk, XList_hk)
  }else{
    if(Method=='iSC-MEB'){
      message("Only use the results in iSC-MEB to remove unwanted variations...")
      HList <- PRECASTObj@resList$iSCMEB$batchEmbed
    }else if(Method=='HarmonyLouvain'){
      message("Only use  the results in Harmony and Louvain to remove unwanted variations...")
      HList <- PRECASTObj@resList$Harmony$batchEmbed
    }else{
      stop("IntegrateSRTData: it does not support this Method!")
    }
    
    
  }
  
  if(!is.null(covariates_use)){
    covariateList <- lapply(PRECASTObj@seulist, function(x) x@meta.data[covariates_use])
  }else{
    covariateList <- NULL
  }
  
  if(is.null(seuList_raw) && !is.null(PRECASTObj@seuList)){
    seuList_raw <- PRECASTObj@seuList
  }
  
  
  
  if(is.null(seuList_raw)){ ## remove the unwanted variation for the selected genes
    defAssay_vec <- sapply(PRECASTObj@seulist, DefaultAssay)
    if(any(defAssay_vec!=defAssay_vec[1])) warning("IntegrateSpaData: there are different default assays in PRECASTObj@seulist that will be used to integrating!")
    n_r <- length(defAssay_vec)
    XList <- lapply(1:n_r,  function(r) Matrix::t(PRECASTObj@seulist[[r]][[defAssay_vec[r]]]@data))
    XList <- lapply(XList, function(x) as.matrix(x))
  }else{ ## remove the unwanted variation for all genes in the data
    ## use the same genes
    
    gene_intersect <- Reduce(intersect, lapply(seuList_raw, row.names))
    ## keep the same spots as the filtered data in PRECASTObj
    seuList_raw <- lapply(1:M, function(r){
      seu_tmp <- seuList_raw[[r]][gene_intersect, ]
      seu_tmp <- seu_tmp[, colnames(PRECASTObj@seulist[[r]])]
      return(seu_tmp)
    })
    defAssay_vec <- sapply(seuList_raw, DefaultAssay)
    if(any(defAssay_vec!=defAssay_vec[1])) warning("IntegrateSpaData: there are different default assays in PRECASTObj@seulist that will be used to integrating!")
    n_r <- length(defAssay_vec)
    XList <- lapply(1:n_r,  function(r) Matrix::t(seuList_raw[[r]][[defAssay_vec[r]]]@data))
    XList <- lapply(XList, function(x) as.matrix(x))
  }
  
  
  .logDiffTime(sprintf(paste0("%s Finish PCA"), "*****"), t1 = tstart, verbose = verbose)
  
  
  
  M <- length(PRECASTObj@seulist)
  p <- ncol(XList[[1]])
  if(is.null(Tm)){ # the temporal information
    Tm <- matrix(0, M, 1)
  }
  
  
  
  
  if(verbose)
    message( "******","Remove the unwanted variations in gene expressions using spatial linear regression...")
  tstart <- Sys.time()
  if(Method=='iSC-MEB'){
    RList <- PRECASTObj@resList$iSCMEB$RList
  }else if(Method=='HarmonyLouvain'){
    clusters_lovain <- PRECASTObj@resList$Louvain$cluster
    cluster_idMat <- model.matrix(~unlist(clusters_lovain)-1)
    RList <- mat2list(cluster_idMat, nvec = nvec)
  }else{
    stop("IntegrateSRTData: it does not support this Method!")
  }
  
  tic <- proc.time()
  correct_List_pois <- correct_genesR(XList, RList, HList, Tm, PRECASTObj@AdjList,
                                      covariateList=covariateList, verbose=verbose)
  toc <- proc.time()
  time_all_pois <- toc[3] - tic[3]
  .logDiffTime(sprintf(paste0("%s Finish unwanted variation removal"), "*****"), t1 = tstart, verbose = verbose)
  
  
  if(verbose)
    message( "******","Sort the results into a Seurat object...")
  tstart <- Sys.time()
  
  XList_correct_all_pois <- correct_List_pois$XList_correct
  ## Add row names and colnames (symbol) for each slice
  for(m in 1:M){
    row.names( XList_correct_all_pois[[m]]) <- paste0("slice", m, "_", row.names( XList[[m]]) )
    colnames(XList_correct_all_pois[[m]]) <- colnames(XList[[m]])
  }
  hX_pois <- matlist2mat(XList_correct_all_pois)
  
  meta_data <- data.frame(batch=factor(get_sampleID(XList)))
  row.names(meta_data) <- row.names(hX_pois)
  rm(XList)
  
  count <- sparseMatrix(i=1,j=1, x=1, dims=dim(t(hX_pois)))
  row.names(count) <- colnames(hX_pois)
  colnames(count) <- row.names(hX_pois)
  seuInt <- CreateSeuratObject(counts = count, project = "ProFAST", assay = "ProFAST", meta.data = meta_data)
  row.names(hX_pois) <- colnames(seuInt)
  seuInt[["ProFAST"]]@data <- t(hX_pois)
  
  seuInt <- Add_embed(matlist2mat(PRECASTObj@resList$ProFAST$hV), seuInt, embed_name = 'profast', assay='ProFAST')
  if(Method=="iSC-MEB"){
    seuInt$cluster <- factor(unlist(PRECASTObj@resList$iSCMEB$cluster))
    seuInt <- Add_embed(matlist2mat(PRECASTObj@resList$iSCMEB$alignedEmbed), seuInt, embed_name = 'iscmeb', assay='ProFAST')
    Idents(seuInt) <- factor(seuInt$cluster)
  }else if(Method=='HarmonyLouvain'){
    seuInt$cluster <-factor(as.numeric(unlist(PRECASTObj@resList$Louvain$cluster)))
    seuInt <- Add_embed(matlist2mat(PRECASTObj@resList$Harmony$harmonyembed), seuInt, embed_name = 'harmony', assay='ProFAST')
    Idents(seuInt) <- factor(as.numeric(unlist(PRECASTObj@resList$Louvain$cluster)))
  }
  
  posList <- lapply(PRECASTObj@seulist, function(x) cbind(x$row, x$col))
  seuInt <- Add_embed(matlist2mat(posList), seuInt, embed_name = 'position', assay='ProFAST')
  
  .logDiffTime(sprintf(paste0("%s Finish results arrangement"), "*****"), t1 = tstart, verbose = verbose)
  
  
  return(seuInt)
}
