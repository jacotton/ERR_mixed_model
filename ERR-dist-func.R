#Breaking the margeff scrip to return posterior distribution


ERRdist <- function(modlist, data, nm="treatment", marg="UID", thresh=0.90)
{
  #Number of chains run
  nchain <- length(modlist)
  #Number of fixed effect parameters
  nfe <- modlist[[1]]$Fixed$nfl[1]
  #Unclear- NULL value
  beta <- b <- vector("list", nchain)
  #Effects in model (fixed and random)
  nms <- colnames(modlist[[1]]$Sol)
  # extracts all RE coeffcients of treatment at any hierarchical level
  cols <- nfe + grep(nm, nms[-c(1:nfe)]) 
  
  # extract coefficients of FEs and REs
  for (i in 1:nchain)
  {
    beta[[i]] <- as.matrix(modlist[[i]]$Sol[, 1:nfe])
    b[[i]] <- as.matrix(modlist[[i]]$Sol[, cols])
  }
  
  # extract design matricies and bind coefficients
  #Matrix of fixed effect values for inidividuals
  X <- modlist[[1]]$X
  #Matrix of random effects values for counts
  Z <- modlist[[1]]$Z
  beta <- do.call(rbind, beta)
  b <- do.call(rbind, b)
  
  # extract fixed effects design matrix columns indicative of treatment
  X <- X[,grep(nm, colnames(X))] 
  # extract random effects design matrix columns indicative of treatment
  Z <- Z[,grep(nm, colnames(Z))]
  # extract coefficients indicative of treatment
  beta <- beta[, grep(nm, colnames(beta))]
  
  # indicate unique rows indicative of single individuals (seems to be +1 at start)
  r <- which( duplicated(cbind(as.matrix(X), as.matrix(Z))) == F )[-1]
  
  # extract relevant rows of design matricies
  Xr <- as.matrix(X[r,])
  Zr <- as.matrix(Z[r,])
  
  # extract indicator of chosen variable (=marg) to marginalize over
  MV <- get(marg, data)[r]
  
  # number of iterations to define length of storage lists
  nit <- nrow(beta)
  logRR <- vector("list", nit)
  ERR <- vector("list", nit)
  
  # calculate ERRs for each iteration
  for (i in 1:nit)
  {
    logRR[[i]] <- (beta[i,]%*%t(X)+b[i,]%*%t(Z))[r] 
  }
  
  # bind togeter ERRs
  logRR <- do.call(cbind, logRR)
  
  # split logRR into mini data frames defined by chosen MV
  logRRsplt <- split( as.data.frame(logRR), MV )
  
  # calculate marginal posteriors distributions of ERR
  # within each group defined by chosen MV 
  # and then report summary statistics of this posterior
  len <- length(logRRsplt)
  newnm <- names(logRRsplt)
  poststats <- vector("list", len)
  
  for (i in 1:len)
  {
    # force ERR to be an array for the case of a single row per individual
    tmpdims <- dim(logRRsplt[[i]])
    ERR <- array( apply(logRRsplt[[i]], 2, FUN = function(x) 1 - exp(x)), dim = tmpdims )
    # coverts each evaluation at posterior to a binary variable
    # above or below chose threshold
    ERRbin <- array( apply(ERR,2, function(x) x<thresh), dim = tmpdims) 
    
    ################# posterior distributions #####################
    ERRpost <- apply(ERR,2, mean)
    ERRthresh <- apply(ERRbin,2, mean)
    
    ################## summary statistics #######################
    poststats[[i]] <- (ERRpost= ERRpost)
  }
  
  data.frame(do.call(rbind, poststats), id=newnm)
  
}
