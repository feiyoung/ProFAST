# # # Model implementation ----------------------------------------------------
# # q = 10
# # sigmaW= c(0.5,0.8,1);
# # sigmaZ = c(1,2, 0.5);
# # qvec = rep(2, 3); # dimension of lantent features of batch effects
# # require(MASS)
# # widthvec <- c(15, 10, 10)
# # heightvec <- c(15, 10, 10)
# # n_vec <- widthvec * heightvec
# #
# # beta_vec <- c(0.8, 1.0, 1.2)
# #
# # K <- 7
# # ### generate class label y from potts model for each sample
# # library(GiRaF)
# # n_vec2 <- rep(n_vec, each=2)
# # heightvec2 <- rep(heightvec, each=2)
# # widthvec2 <- rep(widthvec, each=2)
# # beta_vec2 <- rep(beta_vec, each=2)
# # yList <- list()
# # for(r in 1: length(n_vec2)){
# #   message("r = ", r)
# #   set.seed(r)
# #   y <- sampler.mrf(iter = n_vec2[r], sampler = "Gibbs", h = heightvec2[r], w = widthvec2[r],
# #                    ncolors = K, nei = 4, param = beta_vec2[r],
# #                    initialise = FALSE, view = F)
# #   yList[[r]] <- c(y) + 1
# # }
# # M <- length(n_vec2)
# # sapply(yList, table)
# # ## get the spatial coordinates
# # posList <- lapply(1:M, function(r){
# #   cbind(rep(1:heightvec2[r], widthvec2[r]), rep(1:heightvec2[r], each=widthvec2[r]))
# # })
# # getAdj_reg <- DR.SC:::getAdj_reg
# # ## evaluate the adjecence matrix
# # AdjList <- lapply(posList, getAdj_reg, platform = 'ST')
# #
# #
# # ### generate latent feature V from CAR model for each sample
# # cor.mat <- function(p, rho, type='toeplitz'){
# #
# #   mat <- diag(p)
# #   if(type=='toeplitz'){
# #     for(i in 2:p){
# #       for(j in 1:i){
# #         mat[i,j] <- mat[j,i] <- rho^(abs(i-j))
# #       }
# #     }
# #   }
# #   if(type=='identity'){
# #     mat[mat==0] <- rho
# #   }
# #   return(mat)
# # }
# # cov.mat <- function(sdvec, rho, type='toeplitz'){
# #   p <- length(sdvec)
# #   cormat <- cor.mat(p, rho, type)
# #   covmat <- matrix(0, p, p)
# #   for(i in 1:p){
# #     for(j in 1:p){
# #       covmat[i,j] <- sdvec[i]*sdvec[j]*cormat[i,j]
# #     }
# #   }
# #   return(covmat)
# # }
# # library(LaplacesDemon)
# # library(Matrix)
# #
# # normalizeMatrix <- function(adj){
# #   require(Matrix)
# #   s <- colSums(adj)
# #   n <- nrow(adj)
# #   adj_norm <- adj
# #   for(i in 1:n){
# #     adj_norm[i, i:n] <- adj[i, i:n]/ (2*s[i])
# #     adj_norm[i:n, i] <-  adj_norm[i, i:n]
# #   }
# #   return(adj_norm)
# # }
# # meanMat <- rbind(rep(0, q), rep(c(0,0), length=q), rep(c(0, 0), length=q))
# # meanMat <- rbind(meanMat, meanMat)
# #
# # VList <- list()
# # for(r in 1: 6){
# #   # r <- 1
# #   cat("r = ", r, "\n")
# #   set.seed(r)
# #   M <- matrix(meanMat[r, ], nrow=n_vec2[r], ncol=q, byrow=T)
# #   W <- normalizeMatrix(AdjList[[r]])
# #   U <- as(diag(rep(1, n_vec2[r])), "sparseMatrix") -  W
# #   U2 <- qr.solve(U)
# #   V <- cov.mat(rep(1, q), rho=2*0.2)
# #   U1 <- round(as.matrix(U2),4)
# #   VList[[r]] <- rmatrixnorm(M, U1, V)
# # }
# #
# # colMeans(VList[[3]] )
# # # ### mean-zero autoregressive model
# # # Umat <- matrix(NA, 3, q)
# # # Umat[1,] <- rnorm(q)
# # # Umat[2,] <- 0.5 * Umat[1,] + rnorm(q, sd=0.1)
# # # Umat[3,] <- 0.5 * Umat[2,] + rnorm(q, sd=0.2)
# #
# #
# #
# # p <- 50
# # q = 10
# # sigmaW=c(0.5,0.8,1);
# # sigmaZ = c(1,2, 0.5);
# # qvec=rep(2, 3); # dimension of lantent features of batch effects
# # require(MASS)
# #
# #
# #
# # ## generate deterministic parameters, fixed after generation
# # set.seed(5)
# # sigma2 <- 1
# # LamMat <- rbind(sigma2*(1 + 1.5 * abs(rnorm(p, sd=3))),
# #                 sigma2*(1 + 1*(runif(p))),
# #                 sigma2*(1 + 2*(runif(p))))
# #
# #
# #
# # W <- matrix(rnorm(p*q), p, q)
# # W <- qr.Q(qr(W))
# #
# #
# # M <- 6
# # tvec <- rep(1:3, each=2)
# # ## heter covariance components
# # diagmat = array(0, dim = c(q, q, K))
# # for(k in 1:K){
# #   diag(diagmat[,,k]) <- 1
# # }
# # diag(diagmat[,,1]) = c(10,rep(1,q-1))
# # diag(diagmat[,,2]) = c(1,10,rep(1,q-2))
# # diag(diagmat[,,3]) = c(rep(1,2),10, rep(1,q-3))
# #
# #
# #
# #
# # Mu <- 2*rbind(rep(3, q), rep(0, q), rep(-3, q), rep(2, q), rep(1, q), rep(-1, q), rep(-2, q))
# # Sigma <- diagmat
# #
# # i <- 30
# # set.seed(i)
# #
# # ## generate low-dimensional embedding with biological effects
# # BatchList <- list()
# # q_batch <- 2
# # Zlist <- list()
# # for(r in 1:M){
# #   Z_tmp <- matrix(0, n_vec2[r], q)
# #   for(k in 1:K){
# #     nk <- sum(yList[[r]]==k)
# #     if(nk > 0)
# #       Z_tmp[yList[[r]]==k, ] <- MASS::mvrnorm(nk, Mu[k,],Sigma[,,k])
# #   }
# #   BatchList[[r]] <- MASS::mvrnorm(n_vec2[r], rep(2*r, q_batch), r* diag(rep(1, q_batch)))
# #   Zlist[[r]] <- Z_tmp
# # }
# #
# # # alpha0 <- rep(0, length=p)
# # # alpha0 <- alpha0 - W%*%qr.solve(t(W)%*% W) %*% t(W) %*% alpha0 ## To make alpha0 identifiable
# # beta0 <- rep(0.5, length=q)
# #
# # LamMat2 <- rbind(LamMat, LamMat+ matrix(runif(3*p, 0, 0.2), 3, p))
# # VListall <- list()
# # WmList <- list()
# # # generate log counts
# # XList <- list()
# # for(r in 1: 6){
# #
# #   message("r = ", r)
# #   VListall[[r]] <- Zlist[[r]] + VList[[r]]
# #   tmpMat <- matrix(rnorm(p* q_batch), p, q_batch)
# #   WmList[[r]] <- qr.Q(qr(lm(tmpMat ~ W)$residuals))
# #
# #   XList[[r]] <- (matrix(beta0*tvec[r], nrow=n_vec2[r], ncol=q ,byrow = T) +  VListall[[r]]) %*% t(W) +
# #     + BatchList[[r]] %*% t(WmList[[r]])+ MASS::mvrnorm(n_vec2[r], rep(0,p), diag(LamMat2[r,]))
# #
# # }
# #
# # str(XList)
# #
# #
# #
# #
# #
# #
# #
# # # Test makerDRremovebatch.cpp ---------------------------------------------
# # library(ProFAST)
# # XList_copy <- XList
# # reslistV_T <- ProFAST_run(XList, AdjList = AdjList, q=10, maxIter = 25)
# #
# # str(reslistV_T)
# # head(reslistV_T$hV)
# # XList_count <- lapply(XList, abs)
# # str(XList_count)
# # reslist_count_T <- ProFAST_run(XList_count, AdjList = AdjList,
# #                                           q=10, fit.model = "poisson",maxIter = 25, Vint_zero = T)
# # str(reslist_count_T)
# #
# #
# # head(reslist_count$Mu)
# #
# # cancor(matlist2mat(reslistV$hV), matlist2mat(reslist_count$hV))$cor
# #
# # resList <- fit.iscmeb(reslistV$hV, AdjList, K=5:6, maxIter=25, dr.method = "SpatialFactor",
# #                       coreNum = 2)
# # str(resList)
# # str(resList@paramList)
# # reslist <- selectmodel_iscmeb(resList, c_penalty = 0.3, K=5:6)
# # str(reslist)
# #
# #
# # resList1 <- selectmodel_iscmeb(resList, c_penalty = 0.3, K=5:6, return_fitList = T)
# # str(resList1@paramList)
# 
# 
# 
# # # Test high-level functions using PRECASTObject ---------------------------
# #
# # library(PRECAST)
# # example(CreatePRECASTObject)
# # PRECASTObj
# # PRECASTObj <- AddAdjList(PRECASTObj, platform = "ST")
# # PRECASTObj@AdjList
# # PRECASTObj@parameterList
# #
# # PRECASTObj2 <- AddParSettingProFAST(PRECASTObj)
# #
# # PRECASTObj <- ProFAST(PRECASTObj = PRECASTObj2, fit.model="gaussian")
# # str(PRECASTObj@resList$ProFAST)
# # ## select the number of clusters
# # PRECASTObj <- SelectclustNumber(PRECASTObj)
# #
# # ## Run iSCMEB
# # PRECASTObj <- RuniSCMEB(PRECASTObj)
# # str(PRECASTObj@resList$iSCMEB)
# #
# #
# #
# #
# # # Integration anlaysis
# # seuInt <- IntegrateSRTData(PRECASTObj, seulist_HK=PRECASTObj@seulist, seuList_raw=NULL, covariates_use=NULL, verbose=TRUE)
# #
# # seuInt <- IntegrateSRTData(PRECASTObj, seulist_HK=PRECASTObj@seulist, seuList_raw=NULL, covariates_use="row", verbose=TRUE)
# # library(Seurat)
# # seuInt <- RunTSNE(seuInt, reduction = "alignedprofast")
# # DimPlot(seuInt)
# # seuInt
# #
# #
# 
# 
# 
# # Using simulated data ----------------------------------------------------
# CAR_model_generator <- function(adj, n, r, rho = 0.2, seed.num = 1, phi = 0.99) {
#   set.seed(as.integer(seed.num))
#   message("generating CAR model.\n")
#   M    <- matrix(0, nrow = n, ncol = r)
#   s    <- as(diag(Matrix::colSums(adj)), "sparseMatrix")
#   U    <- qr.solve(s - adj*phi)
#   U    <- ( U + t(U) )/2 # the output of qr.solve fall to through is.positive.definite test
#   Psi  <- outer(1:r, 1:r, FUN = function(x, y) rho^(abs(x - y)))
#   V    <- LaplacesDemon::rmatrixnorm(M, U, Psi)
#   return(
#     list(
#       V = V,
#       Psi = Psi
#     )
#   )
# }
# 
# generate_param <- function(seed = 1, p=2000) {
#   require(MASS)
# 
#   #p <- 2000
#   q <- 10
#   # sigmaW <- c(0.5, 0.8, 1)
#   sigmaZ <- 2 * c(1,2, 0.5)
#   qvec <- rep(2, 3) # dimension of lantent features of batch effects
#   # widthvec <- c(55, 60, 50)
#   # heightvec <- c(55, 50, 60)
#   widthvec <- c(25, 20, 30)
#   heightvec <- c(25, 20, 30)
#   n_vec <- widthvec * heightvec
#   K <- 7
#   mu_value <- 10
#   beta <- c(0.8, 1.0, 1.2)
#   sigma2 <- 2
#   phi <- 0.99
#   rho <- 0.2
# 
#   set.seed(seed)
#   yList <- lapply(as.list(1:3), function(i) GiRaF::sampler.mrf(
#     iter = n_vec[i], sampler = "Gibbs", h = heightvec[i], w = widthvec[i], ncolors = K,
#     nei = 4, param = beta[i], initialise = FALSE, view = FALSE))
#   yList <- lapply(as.list(1:3), function(i) as.vector(yList[[i]]) + 1)
#   posList <- lapply(as.list(1:3), function(i) {
#     cbind(rep(1:heightvec[i], widthvec[i]), rep(1:heightvec[i], each=widthvec[i]))
#   })
#   adjList <- lapply(as.list(1:3), function(i) PRECAST::getAdj_reg(posList[[i]], "ST"))
#   VList <- lapply(as.list(1:3), function(i) {
#     CAR_model_generator(adjList[[i]], n_vec[i], q, rho, i+20000, phi)$V
#   })
# 
#   set.seed(seed)
# 
#   LamMat <- rbind(
#     sigma2*(1 + 1*(runif(p))),
#     sigma2*(1 + 1*(runif(p))),
#     sigma2*(1 + 1*(runif(p))))
# 
#   # generate W (p by q) and Wlist (p by qvec)
#   W <- matrix(rnorm(p * q), p, q)
#   W <- qr.Q(qr(W))
# 
#   # generate mu
#   mu <- matrix(c(rep(0,q),
#                  rep(0,q),
#                  rep(0,q),
#                  rep(0,q),
#                  rep(0,q),
#                  rep(0,q),
#                  rep(0,q)), ncol = K)
#   mu[1,1] = mu_value
#   mu[2,2] = mu_value
#   mu[3,3] = mu_value
#   mu[4,4] = mu_value
#   mu[5,5] = -mu_value
#   mu[6,6] = -mu_value
#   mu[7:q,7] = -mu_value
# 
#   ## heter covariance components
#   diagmat = array(0, dim = c(q, q, K))
#   for(k in 1:K){
#     diag(diagmat[,,k]) <- 1
#   }
#   diag(diagmat[,,1]) = c(10,rep(1,q-1))
#   diag(diagmat[,,2]) = c(1,10,rep(1,q-2))
#   diag(diagmat[,,3]) = c(rep(1,2),10, rep(1,q-3))
#   diag(diagmat[,,4]) = c(rep(1,3),10, rep(1,q-4))
# 
#   Mu <- t(mu)
#   Sigma <- diagmat
#   tauMat <- matrix(0, 3, q)
#   tauMat[2,1] <- 10; tauMat[2,2] <- -10
#   tauMat[3, ] <- rep(5, q);
# 
#   set.seed(seed)
#   diagmat = array(0, dim = c(q, q, 3))
#   for(r in 1:3) {
#     diag(diagmat[,,r]) = abs(rnorm(q, sd = 2))
#   }
#   Psi = diagmat
# 
#   tau0Mat <- matrix(NA, 3, p)
#   for(r in 1:3){
#     set.seed(r+5)
#     tau0 <- rnorm(p, sd=2)
#     tau0Mat[r, ] <- tau0
#   }
# 
#   param <- list(
#     n_vec = n_vec,
#     p = p,
#     q = q,
#     K = K,
#     yList = yList,
#     Mu = Mu,
#     tauMat = tauMat,
#     Psi = Psi,
#     Sigma = Sigma,
#     VList = VList,
#     W = W,
#     LamMat = LamMat,
#     tau0Mat = tau0Mat,
#     posList = posList
#   )
#   return(param)
# }
# 
# generate_seuList_pos <- function(param, sd_signal, sd_batch, seed = 1) {
#   require(MASS)
#   library(Seurat)
#   p <- param$p
#   q <- param$q
#   K <- param$K
#   n_vec <- param$n_vec
#   set.seed(seed)
# 
#   Zlist <- vector("list", 3)
#   for (r in 1:3) {
#     Z_tmp <- matrix(0, n_vec[r], q)
#     for (k in 1:K) {
#       nk <- sum(param$yList[[r]]==k) # include conditional and sequencing batch
#       if (nk > 0) {
#         Z_tmp[param$yList[[r]]==k, ] <- MASS::mvrnorm(
#           nk, param$Mu[k,] + param$tauMat[r,],
#           # param$Sigma[,,k]*sd_signal+sd_batch*(r-1)*diag(q))
#           param$Sigma[,,k]*sd_signal+sd_batch*param$Psi[,,r])
#       }
#     }
#     Zlist[[r]] <- Z_tmp
#   }
#   sapply(Zlist, dim)
# 
#   # generate log counts
#   XtList <- vector("list", 3)
#   for(r in 1:3){
#     X1 <- (Zlist[[r]] + param$VList[[r]] ) %*% t(param$W) +
#       # Zrlist[[r]] %*% t(Wlist[[r]]) +
#       MASS::mvrnorm(n_vec[r], rep(0,p), diag(param$LamMat[r,]))
# 
#     tauMat0 <- matrix(param$tau0Mat[r, ], n_vec[r], p, byrow = T)
#     Eta <- exp((X1 + tauMat0))
#     summary(colSums(Eta))
#     #X1 <- matrix(rpois(n_vec[r]*p, Eta), n_vec[r], p)
#     XtList[[r]] <- matrix(rpois(n_vec[r]*p, Eta), n_vec[r], p)
#   }
# 
#   # create seuList
#   seuList <- vector("list", 3)
#   for (r in 1:3) {
#     count = t(XtList[[r]])
#     y_r = as.factor(param$yList[[r]])
#     levels(y_r) <- c(paste0("Layer",1:6), "WM")
#     meta_data = data.frame(
#       row = param$posList[[r]][,1],
#       col = param$posList[[r]][,2],
#       Group = y_r
#     )
#     rownames(count) = paste0("gene",1:p)
#     colnames(count) = paste0("spot",1:n_vec[r])
#     rownames(meta_data) = colnames(count)
#     seuList[[r]] <- Seurat::CreateSeuratObject(counts=count, meta.data = meta_data)
#   }
# 
#   seuList
# }
# paramList <- generate_param(p=200)
# seuList <- generate_seuList_pos(param = paramList, sd_signal = 2, sd_batch = 0.5)
# seuList
# 
# library(iSC.MEB)
# iSCMEBObj <- CreateiSCMEBObject(seuList = seuList, verbose = FALSE, premin.spots = 0, postmin.spots = 0)
# # simu3 <- seuList
# # save(simu3, file='simu3.rds')
# iSCMEBObj <- CreateNeighbors(iSCMEBObj, platform = "ST")
# ## run PCA to get low dimensional embeddings
# iSCMEBObj <- runPCA(iSCMEBObj, npcs = 15, pca.method = "APCA")
# ## Add a model setting in advance for an iSCMEBObj object. verbose = TRUE helps outputing the
# ## information in the algorithm.
# iSCMEBObj <- SetModelParameters(iSCMEBObj, verbose = TRUE)
# iSCMEBObj <- iSCMEB(iSCMEBObj, K = 7)
# iSCMEBObj <- SelectModel(iSCMEBObj)
# LabelList <- lapply(iSCMEBObj@seulist, function(seu) seu@meta.data$Group)
# ARI <- function(x, y) mclust::adjustedRandIndex(x, y)
# ari_sections <- sapply(1:3, function(i) ARI(idents(iSCMEBObj)[[i]], LabelList[[i]]))
# ari_all <- ARI(unlist(idents(iSCMEBObj)), unlist(LabelList))
# print(ari_sections)
# 
# seuList <- simu3
# library(ProFAST)
# row.names(seuList[[1]])[1:10]
# 
# ## Get the gene-by-spot read count matrices
# countList <- lapply(seuList, function(x) x[["RNA"]]@counts)
# M <- length(countList)
# 
# 
# 
# 
# ## Get the meta data of each spot for each data batch
# metadataList <- lapply(seuList, function(x) x@meta.data)
# 
# for (r in 1:M) {
#   meta_data <- metadataList[[r]]
#   all(c("row", "col") %in% colnames(meta_data))  ## the names are correct!
#   print(head(meta_data[, c("row", "col")]))
# }
# ## ensure the row.names of metadata in metaList are the same as that of colnames count matrix
# ## in countList
# 
# for (r in 1:M) {
#   row.names(metadataList[[r]]) <- colnames(countList[[r]])
# }
# 
# 
# ## Create the Seurat list object
# 
# simu2 <- list()
# for (r in 1:M) {
#   simu2[[r]] <- CreateSeuratObject(counts = countList[[r]], meta.data = metadataList[[r]], project = "SimuProFAST")
# }
# 
# library(PRECAST)
# library(ProFAST)
# ## Create PRECASTObject.
# set.seed(2023)
# PRECASTObj <- CreatePRECASTObject(simu2, project = "simu2", gene.number =180, selectGenesMethod = "HVGs",
#                                   premin.spots = 20, premin.features = 20, postmin.spots = 1, postmin.features = 1)
# 
# #### Run ProFAST
# PRECASTObj <- AddAdjList(PRECASTObj, platform = "ST")
# 
# ### Add basic settings for ProFAST
# PRECASTObj <- AddParSettingProFAST(PRECASTObj)
# unlist(PRECASTObj@parameterList) # check the parameter settings
# 
# PRECASTObj <- ProFAST(PRECASTObj = PRECASTObj, fit.model="poisson")
# #PRECASTObj@seulist[[1]]$Group
# yList <- lapply(PRECASTObj@seulist, function(x) x$Group)
# ### Evaluate the MacR2
# (MacVec_gauss <- evaluate_DR_PF2(PRECASTObj@resList$ProFAST$hV, yList))
# #(MacVec <- evaluate_DR_PF2(PRECASTObj@resList$ProFAST$hV, yList))
# 
# 
# str(PRECASTObj@resList$ProFAST)
# ## select the number of clusters
# PRECASTObj <- RunHarmonyLouvain(PRECASTObj, resolution = 0.4)
# lapply(1:M, function(r) mclust::adjustedRandIndex(PRECASTObj@resList$Louvain$cluster[[r]], yList[[r]]))
# 
# 
# ## Run iSCMEB
# PRECASTObj <- RuniSCMEB(PRECASTObj, seed=2023)
# str(PRECASTObj@resList$iSCMEB)
# lapply(1:M, function(r) mclust::adjustedRandIndex(PRECASTObj@resList$iSCMEB$cluster[[r]], yList[[r]]))
# 
# 
# 
# # Integration anlaysis
# ## select the HK genes
# HKgenes <- SelectHKgenes(simu2, species= "Human", HK.number=200)
# seulist_HK <- NULL
# seuInt <- IntegrateSRTData(PRECASTObj, seulist_HK=NULL, Method = "HarmonyLouvain", seuList_raw=NULL, covariates_use=NULL, verbose=TRUE)
# 
# head(PRECASTObj@seulist[[1]]@meta.data)
# ### try covariates
# IntegrateSRTData(PRECASTObj, seulist_HK=NULL, Method = "HarmonyLouvain", 
#                  seuList_raw=NULL, covariates_use=c("row", "col"), verbose=TRUE)
# 
# IntegrateSRTData(PRECASTObj, seulist_HK=NULL, Method = "HarmonyLouvain", 
#                  seuList_raw=NULL, covariates_use=c("Group", "col"), verbose=TRUE)
# 
# IntegrateSRTData(PRECASTObj, seulist_HK=NULL, Method = "iSC-MEB", 
#                  seuList_raw=NULL, covariates_use=c("Group", "col"), verbose=TRUE)
# IntegrateSRTData(PRECASTObj, seulist_HK=NULL, Method = "iSC-MEB", 
#                  seuList_raw=seuList, covariates_use=c("row", "col"), verbose=TRUE)
# 
# # seuInt <- IntegrateSRTData(PRECASTObj, seulist_HK=PRECASTObj@seulist, seuList_raw=dlpfc2)
# library(Seurat)
# seuInt <- RunTSNE(seuInt, reduction = "harmony")
# seuInt <- AddTSNE(seuInt, n_comp=3, reduction = 'harmony', assay = 'ProFAST')
# #seuInt <- RunTSNE(seuInt, reduction = "alignedprofast")
# cols_cluster <- chooseColors(palettes_name = "Nature 10", n_colors = 7, plot_colors = TRUE)
# p_batch <- dimPlot(seuInt, item='batch')
# p_cluster <- dimPlot(seuInt, item='cluster', cols = cols_cluster)
# drawFigs(list(p_batch, p_cluster), layout.dim = c(1,2))
# 
# DimPlot(seuInt)
# 
# SpaPlot(seuInt, item = "cluster", batch = NULL, point_size = 1, cols = cols_cluster, combine = T,
#                nrow.legend = 8)
# 
# 
# #iSC-MEB ---------------------------------------
# ## Visualize the batch effects
# seuInt <- AddTSNE(seuInt, n_comp=2, reduction = 'profast', assay = 'ProFAST')
# tsne2_profastp <- seuInt@reductions$tSNE
# p_batch_profastp <- dimPlot(seuInt, item='batch')
# seuInt$cluster_true <- factor(unlist(yList))
# p_cluster_profastp <- dimPlot(seuInt, item='cluster_true', cols = cols_cluster)
# drawFigs(list(p_batch_profastp, p_cluster_profastp), layout.dim = c(1,2))
# 
# 
# ### see Harmony
# seuInt <- AddTSNE(seuInt, n_comp=2, reduction = 'harmony', assay = 'ProFAST')
# tsne2_harmony <- seuInt@reductions$tSNE
# p_batch_harmony <- dimPlot(seuInt, item='batch')
# p_cluster_harmony <- dimPlot(seuInt, item='cluster_true', cols = cols_cluster)
# drawFigs(list(p_batch_harmony, p_cluster_harmony), layout.dim = c(1,2))
# 
# 
# seuInt2 <- IntegrateSRTData(PRECASTObj, seulist_HK=NULL, Method = "iSC-MEB", seuList_raw=NULL, covariates_use=NULL, verbose=TRUE)
# ## visualize the embed of iscmeb
# seuInt2 <- AddTSNE(seuInt2, n_comp=2, reduction = 'iscmeb', assay = 'ProFAST')
# tsne2_iscmeb <- seuInt2@reductions$tSNE
# p_batch_iscmeb <- dimPlot(seuInt2, item='batch')
# seuInt2$cluster_true <- factor(unlist(yList))
# p_cluster_iscmeb <- dimPlot(seuInt2, item='cluster_true', cols = cols_cluster)
# drawFigs(list(p_batch_iscmeb, p_cluster_iscmeb), layout.dim = c(1,2))
# # Using the real DLPFC data -----------------------------------------------
# # dat <- NULL
# # get_r2_mcfadden <- function(embeds, y){
# #   library(nnet)
# #   library(performance)
# #   y <- as.numeric(as.factor(y))
# #   hq <- ncol(embeds)
# #   dat <- as.data.frame(cbind(y=y, x=embeds))
# #   dat$y <- factor(dat$y)
# #   name <-  c('y', paste0('V', 1:hq))
# #   names(dat) <-name
# #   formu <- paste0("y~")
# #   for(i in 1:hq){
# #     if(i < hq){
# #       formu <- paste(formu, name[i+1], seq='+')
# #     }else{
# #       formu <- paste(formu, name[i+1], seq='')
# #     }
# #
# #   }
# #   model1 <- nnet::multinom(as.formula(formu), data = dat)
# #   R2 <- r2_mcfadden(model1)
# #   return(R2$R2_adjusted)
# # }
# evaluate_DR_PF2 <- function(hZList, yList){
#   MacVec_pro_hV <- rep(NA, length(yList))
#   for(r in 1:length(yList)){
#     #r <- 1
#     MacVec_pro_hV[r] <- get_r2_mcfadden(hZList[[r]], yList[[r]])
#   }
# 
#   MacVec_pro_hV
# }
# 
# seuList <- readRDS("./vignettes/seulist2_ID9_10.RDS")
# 
# row.names(seuList[[1]])[1:10]
# 
# ## Get the gene-by-spot read count matrices
# countList <- lapply(seuList, function(x) x[["RNA"]]@counts)
# M <- length(countList)
# ### transfer the Ensembl ID to symbol
# for(r in 1:M){
#   row.names(countList[[r]]) <- transferGeneNames(row.names(countList[[r]]), now_name = "ensembl",
#                                                  to_name="symbol",
#                                                  species="Human", Method='eg.db')
# }
# 
# 
# 
# ## Get the meta data of each spot for each data batch
# metadataList <- lapply(seuList, function(x) x@meta.data)
# 
# for (r in 1:M) {
#   meta_data <- metadataList[[r]]
#   all(c("row", "col") %in% colnames(meta_data))  ## the names are correct!
#   print(head(meta_data[, c("row", "col")]))
# }
# 
# ## ensure the row.names of metadata in metaList are the same as that of colnames count matrix
# ## in countList
# 
# for (r in 1:M) {
#   row.names(metadataList[[r]]) <- colnames(countList[[r]])
# }
# 
# 
# ## Create the Seurat list object
# 
# seuList <- list()
# for (r in 1:M) {
#   seuList[[r]] <- CreateSeuratObject(counts = countList[[r]], meta.data = metadataList[[r]], project = "BreastCancerPRECAST")
# }
# 
# dlpfc2 <- seuList
# 
# head(meta_data[, c("row", "col")])
# 
# ## Create PRECASTObject.
# set.seed(2023)
# PRECASTObj <- CreatePRECASTObject(dlpfc2, project = "DLPFC2", gene.number = 2000, selectGenesMethod = "HVGs",
#                                   premin.spots = 20, premin.features = 20, postmin.spots = 1, postmin.features = 1)
# 
# #### Run ProFAST
# PRECASTObj <- AddAdjList(PRECASTObj, platform = "Visium")
# 
# ### Add basic settings for ProFAST
# PRECASTObj <- AddParSettingProFAST(PRECASTObj)
# unlist(PRECASTObj@parameterList) # check the parameter settings
# 
# PRECASTObj <- ProFAST(PRECASTObj = PRECASTObj, fit.model="gaussian")
# 
# yList <- lapply(PRECASTObj@seulist, function(x) x$layer_guess_reordered)
# colMeans(PRECASTObj@resList$ProFAST$hV[[1]])
# ### Evaluate the MacR2
# (MacVec_gauss <- evaluate_DR_PF2(PRECASTObj@resList$ProFAST$hV, yList))
# (MacVec <- evaluate_DR_PF2(PRECASTObj@resList$ProFAST$hV, yList))
# 
# 
# str(PRECASTObj@resList$ProFAST)
# ## select the number of clusters
# PRECASTObj <- RunHarmonyLouvain(PRECASTObj, resolution = 0.3) ## six clusters
# lapply(1:M, function(r) mclust::adjustedRandIndex(PRECASTObj@resList$Louvain$cluster[[r]], yList[[r]]))
# 
# ## Run iSCMEB
# PRECASTObj <- RuniSCMEB(PRECASTObj, seed=1)
# str(PRECASTObj@resList$iSCMEB)
# lapply(1:M, function(r) mclust::adjustedRandIndex(PRECASTObj@resList$iSCMEB$cluster[[r]], yList[[r]]))
# 
# 
# # Integration anlaysis
# ## select the HK genes
# HKgenes <- SelectHKgenes(dlpfc2, species= "Human", HK.number=200)
# seulist_HK <- lapply(dlpfc2, function(x) x[HKgenes, ])
# seuInt <- IntegrateSRTData(PRECASTObj, seulist_HK=seulist_HK,
#                            seuList_raw=NULL, covariates_use=NULL, verbose=TRUE)
# 
# seuInt <- IntegrateSRTData(PRECASTObj, seulist_HK=seulist_HK, Method = "HarmonyLouvain",
#                            seuList_raw=NULL, covariates_use=NULL, verbose=TRUE)
# # seuInt <- IntegrateSRTData(PRECASTObj, seulist_HK=PRECASTObj@seulist, seuList_raw=dlpfc2)
# library(Seurat)
# #seuInt <- RunTSNE(seuInt, reduction = "alignedprofast")
# seuInt <- RunTSNE(seuInt, reduction = "harmony")
# # DimPlot(seuInt)
# # seuInt
# # Idents(seuInt)
# 
# p_batch <- dimPlot(seuInt, item='batch')
# p_cluster <- dimPlot(seuInt, item='cluster')
# drawFigs(list(p_batch, p_cluster), layout.dim = c(1,2))
# 
# cols_cluster <- chooseColors(palettes_name = "Nature 10", n_colors = 6, plot_colors = TRUE)
# p12 <- SpaPlot(seuInt, item = "cluster", batch = NULL, point_size = 1, cols = cols_cluster, combine = FALSE,
#                nrow.legend = 8)
# library(ggplot2)
# p12 <- lapply(p12, function(x) x+ coord_flip() + scale_x_reverse())
# drawFigs(p12, layout.dim = c(1,2), common.legend = T)
# 
# 
# ### Compare with Seurat V3
# # normalize and identify variable features for each dataset independently
# ifnb.list <- lapply(X = PRECASTObj@seulist, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# 
# # select features that are repeatedly variable across datasets for integration
# features <- SelectIntegrationFeatures(object.list = ifnb.list)
# immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# # this command creates an 'integrated' data assay
# immune.combined <- IntegrateData(anchorset = immune.anchors)
# 
# # specify that we will perform downstream analysis on the corrected data note that the
# # original unmodified data still resides in the 'RNA' assay
# DefaultAssay(immune.combined) <- "integrated"
# 
# npcs <- 15 # use 15 PCs for fair comparison
# # Run the standard workflow for visualization and clustering
# immune.combined <- ScaleData(immune.combined, verbose = FALSE)
# immune.combined <- RunPCA(immune.combined, npcs = npcs, verbose = FALSE)
# immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:npcs)
# immune.combined <- FindClusters(immune.combined, resolution = 0.4)
# ## calculate the ARI
# hyList_seurat <- vec2list(immune.combined$seurat_clusters, nvec=sapply(PRECASTObj@seulist, ncol) )
# lapply(1:M, function(r) mclust::adjustedRandIndex(hyList_seurat[[r]], yList[[r]]))
# 
