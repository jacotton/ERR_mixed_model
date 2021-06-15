## function to calculate marginal efficacy over chosen variable (random or fixed effects)

## modlist is a list of fitted models with different starting values from MCMCglmm
## (i.e. different Markov chains)

## data is the data frame fed to MCMCglmm during model fitting

## nm is the variable that indicates whether an observation was made before (0)
## or after treatment (1)

## marg is the name of the variable over which to marginalize (average)

## thresh defines a threshold for evaluating the probability that 
## a marginalized ERR is below this value (deafult sett to 90%)

margeff <- function(modlist, data, nm="treatment", marg="numID", thresh=0.90)
{
  
  nchain <- length(modlist)
  nfe <- modlist[[1]]$Fixed$nfl[1]
  beta <- b <- vector("list", nchain)
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
  X <- modlist[[1]]$X
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
  
  # extract indicator of variable to marginalize over
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
    mn <- mean(ERRpost)
    med <- median(ERRpost)
    lwr <- as.numeric( quantile(ERRpost, probs=c(0.025)) )
    upr <- as.numeric( quantile(ERRpost, probs=c(0.975)) )
    mn_thresh <- mean(ERRthresh)
    med_thresh <- median(ERRthresh)
    lwr_thresh <- as.numeric( quantile(ERRthresh, probs=c(0.025)) )
    upr_thresh <- as.numeric( quantile(ERRthresh, probs=c(0.975)) )
    poststats[[i]] <- c(mn=mn, med = med, lwr=lwr, upr=upr, 
                        mn_thresh = mn_thresh, med_thresh = med_thresh,
                        lwr_thresh = lwr_thresh, upr_thresh = upr_thresh)
  }
  
  data.frame(do.call(rbind, poststats), id=newnm)
  
}
