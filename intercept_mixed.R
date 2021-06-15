## Modified function to determine egg counts pre-treatment from model
#Seems to work

margeff_intercept <- function(modlist, data, nm="treatment", marg="numID")
{
  #Number of chains run
  nchain1 <- length(modlist)
  #Number of fixed effect parameters
  nfe1 <- modlist[[1]]$Fixed$nfl[1]
  #Unclear- NULL value
  beta1 <- b1 <- vector("list", nchain1)
  #Effects in model (fixed and random)
  nms1 <- colnames(modlist[[1]]$Sol)
  # extracts all RE coeffcients of treatment at any hierarchical level
  cols1 <- nfe1 + grep(nm, nms1[-c(1:nfe1)], invert=TRUE) 
  
  # extract coefficients of FEs and REs
  for (i in 1:nchain1)
  {
    beta1[[i]] <- as.matrix(modlist[[i]]$Sol[, 1:nfe1])
    b1[[i]] <- as.matrix(modlist[[i]]$Sol[, cols1])
  }
  
  # extract design matricies and bind coefficients
  #Matrix of fixed effect values for inidividuals
  X1 <- modlist[[1]]$X
  #Matrix of random effects values for counts
  Z1 <- modlist[[1]]$Z
  beta1 <- do.call(rbind, beta1)
  b1 <- do.call(rbind, b1)
  
  # extract fixed effects design matrix columns indicative of INTERCEPT
  X1 <- X1[,grep(nm, colnames(X1),invert=TRUE)] 
  # extract random effects design matrix columns indicative of INTERCEPT
  Z1 <- Z1[,grep(nm, colnames(Z1),invert=TRUE)]
  # extract coefficients indicative of INTERCEPT
  beta1 <- beta1[, grep(nm, colnames(beta1), invert=TRUE)]
  
  # indicate unique rows indicative of single individuals
  r1 <- which( duplicated(cbind(as.matrix(X1), as.matrix(Z1))) == F )
  
  # extract relevant rows of design matricies
  Xr1 <- as.matrix(X1[r1,])
  Zr1 <- as.matrix(Z1[r1,])
  
  # extract indicator of chosen variable (=marg) to marginalize over
  MV1 <- get(marg, data)[r1]
  
  # number of iterations to define length of storage lists
  nit1 <- nrow(beta1)
  logRR1 <- vector("list", nit1)
  ERR1 <- vector("list", nit1)
  
  # calculate ERRs for each iteration
  for (i in 1:nit1)
  {
    logRR1[[i]] <- (beta1[i,]%*%t(X1)+b1[i,]%*%t(Z1))[r1] 
  }
  
  # bind togeter ERRs
  logRR1 <- do.call(cbind, logRR1)
  
  # split logRR into mini data frames defined by chosen MV
  logRRsplt1 <- split( as.data.frame(logRR1), MV1 )
  
  # calculate marginal posteriors distributions of ERR
  # within each group defined by chosen MV 
  # and then report summary statistics of this posterior
  len1 <- length(logRRsplt1)
  newnm1 <- names(logRRsplt1)
  poststats1 <- vector("list", len1)
  
  for (i in 1:len1)
  {
    # force ERR to be an array for the case of a single row per individual
    tmpdims1 <- dim(logRRsplt1[[i]])
    ERR1 <- array( apply(logRRsplt1[[i]], 2, FUN = function(x) exp(x)), dim = tmpdims1 )
    # coverts each evaluation at posterior to a binary variable
    # above or below chose threshold
    #ERRbin <- array( apply(ERR,2, function(x) x<thresh), dim = tmpdims) 
    
    ################# posterior distributions #####################
    ERRpost1 <- apply(ERR1,2, mean)
    #ERRthresh <- apply(ERRbin,2, mean)
    
    ################## summary statistics #######################
    mn1 <- mean(ERRpost1)
    med1 <- median(ERRpost1)
    lwr1 <- as.numeric( quantile(ERRpost1, probs=c(0.025)) )
    upr1 <- as.numeric( quantile(ERRpost1, probs=c(0.975)) )
    #mn_thresh <- mean(ERRthresh)
    #med_thresh <- median(ERRthresh)
    #lwr_thresh <- as.numeric( quantile(ERRthresh, probs=c(0.025)) )
    #upr_thresh <- as.numeric( quantile(ERRthresh, probs=c(0.975)) )
    poststats1[[i]] <- c(mn=mn1, med = med1, lwr=lwr1, upr=upr1)
  }
  
  data.frame(do.call(rbind, poststats1), id=newnm1)
  
}
