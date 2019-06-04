#####################################################################
#Code to parallelze KFE to Woodchester Park clusters
#parallelization implemented with foreach/dopar
#rm(list = ls())

#Comments:
#02/04/19: fixed t0_start bug
#15/04/19: removed iterpc (not necessary anymore)
#22/04/19: upload new version of pairs_calc

#Author:
#Gianluigi Rossi
#####################################################################
#packages
require(readr)
require(igraph)
require(ggplot2)
require(parallel)
require(doParallel)
require(foreach)
require(deSolve)

#####################################################################
#upload code to solve KFE equations, without printing messages
source("190424_KFE_pairs_calc_noPrints.R")

#####################################################################
#FUNCTION INCLUDING MCMC
KFE_cluster_calc <- function(
  #DATA
  pairs.list = NA,
  #MCMC algorithm
  min_like = -500,
  niter = 500,
  #initial parameters, uniform prior 
  #(since the parallelization algorijtm takes one row only at the time, this values are resampled for each pair)
  betaTB_start = runif(1,0,5)/365,
  latperTB_start = runif(1,1,365*5),
  subInf_start = runif(1,0,.5)/365,
  #not estimated params
  maxLifExp = 365*8,
  #sd for normal distrib resampling (kept wide on purpose)
  sd_betaTB = 0.01,
  sd_latperTB = 0.1,
  sd_subInf = 0.01
){
  
  #total number of pairs
  num.pairs <- nrow(pairs.list)
  t0_start <- NA
  
  #store results
  all.betaTB <- numeric(num.pairs)
  all.lateperTB <- numeric(num.pairs)
  all.subInf <- numeric(num.pairs)
  all.loglike <- numeric(num.pairs)
  all.t0 <- numeric(num.pairs)
  
  ################
  for(PAIR in 1:num.pairs){
 # PAIR <- 1
    if(pairs.list$SameAnimal[PAIR] == 0){
    #STARTING MCMC ALGORITHM
    vpois_loglike = numeric(0)
    
    vBetaTB    <- numeric(0)
    vLatperTB  <- numeric(0)
    vSubInf    <- numeric(0)
    vt0        <- numeric(0)
  
    vBetaTB_mcmc    <- numeric(0)
    vLatperTB_mcmc  <- numeric(0)
    vSubInf_mcmc    <- numeric(0)
    vt0_mcmc        <- numeric(0)
    
    ################  
    for(iter in 1:niter){
    #iter <- 1
      betaTB   <- -1
      latperTB <- -1
      subInf   <- -1
      while(betaTB<=0 | betaTB > 5/365) betaTB   <- rnorm(1,betaTB_start, sd = sd_betaTB)/365
      while(latperTB<=0|latperTB>365*5) latperTB <- rnorm(1,latperTB_start,sd = sd_latperTB)
      while(subInf<=0|subInf>.5/365)    subInf   <- rnorm(1,subInf_start,sd = sd_subInf)
    
      sel.Params = c("tau"        = betaTB,
                     "sigma"      = 1/latperTB,
                     "max.life.exp"   = maxLifExp,
                     "sub.inf"    = subInf,
                     "sub.lat"    = 1)
      
      ################
      COMPUT <- KFE_SEI_pairs_calc(pairs.list = pairs.list[PAIR,],
                                    all.Params = sel.Params,
                                    t0.fixed = NA,
                                    integral.time = 12*10)
      
      ################
      pois_loglike  <- log(max(COMPUT$chain.AC.prob))
      vpois_loglike <- append(vpois_loglike,pois_loglike)
      
      #check if infection is possible
      if(pois_loglike > -Inf){
        t0 <- COMPUT$t0[which.max(COMPUT$chain.AC.prob)]
        vt0           <- append(vt0, t0)
      
        vBetaTB   <- append(vBetaTB, betaTB)
        vLatperTB <- append(vLatperTB, latperTB)
        vSubInf   <- append(vSubInf, subInf)

        ################
        
        Rndm <- runif(1,0,1)
        if(iter==1|exp(pois_loglike-min_like)>Rndm|length(which(is.na(vBetaTB_mcmc))) == length(vBetaTB_mcmc)){
          #cat("New values chosen") 
          min_like     <- pois_loglike
          betaTB_start   <- betaTB
          latperTB_start <- latperTB
          subInf_start   <- subInf
          t0_start       <- t0
        }
    
        vBetaTB_mcmc   <- append(vBetaTB_mcmc, betaTB_start)
        vLatperTB_mcmc <- append(vLatperTB_mcmc, latperTB_start)
        vSubInf_mcmc   <- append(vSubInf_mcmc, subInf_start)
        vt0_mcmc       <- append(vt0_mcmc,t0_start)
    
      }else{
        vBetaTB_mcmc   <- append(vBetaTB_mcmc, NA)
        vLatperTB_mcmc <- append(vLatperTB_mcmc, NA)
        vSubInf_mcmc   <- append(vSubInf_mcmc, NA)
        vt0_mcmc       <- append(vt0_mcmc, NA)
      }
      #cat(paste("\r",iter))
    }
    ################
    if(any(vpois_loglike>-Inf)){
      #plot(vpois_loglike, x = 1:niter, type = "l", col = 4)
      #plot(1/vSubInf, x = 1:niter, type = "l", col = 4)
      #points(1/vSubInf_mcmc, x = 1:niter, type = "l", col = 2)
      ################
      all.betaTB[PAIR]    <- vBetaTB_mcmc[which.max(vpois_loglike)]
      all.lateperTB[PAIR] <- vLatperTB_mcmc[which.max(vpois_loglike)]
      all.subInf[PAIR]    <- vSubInf_mcmc[which.max(vpois_loglike)]
      all.loglike[PAIR]   <- max(vpois_loglike)
      all.t0[PAIR]        <- vt0[which.max(vpois_loglike)]
    }else{
      #print("Infection not possible")
      all.betaTB[PAIR]    <- NA
      all.lateperTB[PAIR] <- NA
      all.subInf[PAIR]    <- NA
      all.loglike[PAIR]   <- NA
      all.t0[PAIR]        <- NA
    }
    }else{
      #sequences coming from same individual
      all.betaTB[PAIR]    <- NA
      all.lateperTB[PAIR] <- NA
      all.subInf[PAIR]    <- NA
      all.loglike[PAIR]   <- NA
      all.t0[PAIR]        <- NA
    }
      
    ################
    #print(PAIR)
    #save.image("190319_pairProbab_WPcS_part.Rdata")
  }
  
  return(list(
    all.betaTB = all.betaTB,
    all.lateperTB = all.lateperTB,
    all.subInf = all.subInf,
    all.loglike = all.loglike,
    all.t0 = all.t0)
    )
}
#test
#prova1 <- KFE_cluster_calc(pairs.list=pairs.list.cl0[1:5,], niter = 10)
#t(matrix(unlist(prova1), nrow =5))

#####################################################################
#PARALLEL COMPUTING
KFE_cluster_calc_par <- function(pairs.List, Cores, Niter, ...){
  
  NumPa <- nrow(pairs.List)
  #print(NumPa)
  registerDoParallel(cores = Cores)
  
  foreach(DL = 1:NumPa) %dopar% {
    #print(DL)
    KFE_cluster_calc(pairs.list=pairs.List[DL,], niter = Niter)
  }
}
#test
#GGG <- KFE_cluster_calc_par(pairs.List = pairs.list.cl0[1:8,], Cores=4, Niter = 5)
#RESU <- data.frame(t(matrix(unlist(GGG),nrow = 6)))
#names(RESU) <- names(unlist(GGG))[1:6]

#####################################################################


