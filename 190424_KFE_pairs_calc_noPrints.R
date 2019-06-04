#################################################
#Function to calculate the KFE for triads data

#Comments:
#Added start time and adjusted to permit 1 time estimations
#12/03/19: removed prints
#05/04/19 -> fixed maximum substitutions bug (targeting right column now)
#14/04/19 -> changed KFEmatrixSetup source
#         -> correct error: if tCp < tAp I move the computation until A is sampled
#         -> fixed the absorbing SNPs status bug
#24/04/19 -> life.exp changed in max.life.exp

#Author:
#Gianluigi Rossi
#################################################
#KFE solving functions
source("190510_KFEmatrixSetup_v16_noPrints.R")

#################################################
#function to solve the matrix
solveKFE <- function(time, init.cond, parms)
  with(parms, {
    P <- init.cond
    dP <- P%*%model.mat
    return(list(c(dP)))
  })

#################################################
KFE_SEI_pairs_calc <- function(
  #table with triad informations
  pairs.list = NA, #
  #parameters (in days)
  all.Params = c("tau" = 1,
                 "sigma" = .5,
                 "max.life.exp" = 365,
                 "sub.inf" = 0.000025,
                 "sub.lat" = 1), #
  #starting time
  t0.fixed = NA,
  #number of initial time steps
  integral.time = 120
){
 
  #timings
  tAp.all <- pairs.list$sample.t1
  tCp.all <- pairs.list$sample.t2
  
  if(is.na(t0.fixed) & integral.time >= 2){
    t0.mat <- apply(matrix(tAp.all, ncol = 1), 1, function(x){seq(from = x - 1 - all.Params["max.life.exp"], to = x -1, by = all.Params["max.life.exp"]/(integral.time-1))})
  }else if(!is.na(t0.fixed) & integral.time >= 2){
    t0.mat <- apply(matrix(tAp.all, ncol = 1), 1, function(x){seq(from = x -1 - t0.fixed, to = x -1, by = (x - 1 - t0.fixed)/(integral.time-1))})
  }else if(integral.time == 1){
    t0.mat <- matrix(t0.fixed, nrow = 1, ncol = nrow(pairs.list))
  }
  #maximum substitutions
  subMax <- max(c(pairs.list$snps.d1, pairs.list$snps.d2)) + 1
  
  #store results all pairs
  p.c2.Tstart.allPairs <- matrix(NA,
                               ncol = ncol(t0.mat),
                               nrow = nrow(t0.mat))

  set.pair2 <- set.KFE.system.SEI(
    n.ind = 2,
    sub.max = subMax,
    params = all.Params,
    simplify.mat = TRUE,
    chain.start = TRUE,
    second.inf = FALSE)
  
  ###
  set.ind.A <- set.KFE.system.SEI(
    n.ind = 1,
    sub.max = subMax,
    params = all.Params,
    simplify.mat = TRUE,
    chain.start = TRUE,
    second.inf = FALSE)
  
  ###
  set.ind.C <- set.KFE.system.SEI(
    n.ind = 1,
    sub.max = subMax,
    params = all.Params,
    simplify.mat = TRUE,
    chain.start = FALSE,
    second.inf = TRUE)
  
  ##
  mat.pair2 <- set.pair2$KFE.matrix
  mat.ind.A <- set.ind.A$KFE.matrix
  mat.ind.C <- set.ind.C$KFE.matrix
  
for(Num in 1:nrow(pairs.list)){ # Num <- 1
  #time A sampling
  tAp <- tAp.all[Num]
  #time C sampling
  tCp <- tCp.all[Num]
  #timeline
  t0.vec <- t0.mat[,Num]
  t0.len <- length(t0.vec)
  
  #store results single pair
  p.c2.Tstart <- numeric(t0.len)
  
  #select divergent SNPs of A
  div.A <- pairs.list$snps.d1[Num]
  div.C <- pairs.list$snps.d2[Num]
  #pairs.list[Num,]
 
  #chain 2 -> A infected C
  initcond2.1 <- numeric(set.pair2$tot.statuses)
  initcond2.1[which(apply(set.pair2$names.statuses.tab[,1:3], 1, function(x) all(x == c("E", 0,0))))] <- 1
  
  #time start
  for(tstart in 1:t0.len){
    #set startng time from vector - tstart <- 88
    t0 <- t0.vec[tstart]

    ###################################    
    if(tAp < tCp){ #A SAMPLED BEFORE C
     
      ##########
      #calculate probabilities until A was sampled at tAp
      CoarseDeno <- 10
      OUT2.P1T <- ode(y = initcond2.1, times = seq(t0,tAp,by=(tAp-t0)/CoarseDeno), func = solveKFE, parms = list(model.mat=mat.pair2), method = "radau")
      lastrow2.1 <- nrow(OUT2.P1T)
      #check for vals < 0
      while(any(OUT2.P1T[lastrow2.1,-1] < 0)){
        if(any(abs(OUT2.P1T[lastrow2.1,which(OUT2.P1T[lastrow2.1,-1] < 0)+1]) < 10^-30) | CoarseDeno > 2000){
          OUT2.P1T[lastrow2.1,which(OUT2.P1T[lastrow2.1,-1] < 0)+1] <- 0
        }else{
          CoarseDeno <- CoarseDeno*2
          OUT2.P1T <- ode(y = initcond2.1, times = seq(t0,tAp,by=(tAp-t0)/CoarseDeno), func = solveKFE, parms = list(model.mat=mat.pair2), method = "radau")
          lastrow2.1 <- nrow(OUT2.P1T)
          }
      }

      ##########
      #chain 2 -> A infected C
      OUT2.P1 <- OUT2.P1T[,which(set.pair2$names.statuses.tab[,3] == div.A & 
                                 set.pair2$names.statuses.tab[,1] == "I") + 1]
      
      names.statuses.tab.pair2.select <- set.pair2$names.statuses.tab[which(set.pair2$names.statuses.tab[,3] == div.A & 
                                                                            set.pair2$names.statuses.tab[,1] == "I"),]
      
      #set initial conditions for second step
      initcond2.2 <- numeric(set.ind.C$tot.statuses)
      for(stAt2 in 1:set.ind.C$tot.statuses){
        pos.to.calc2 <- which(apply(names.statuses.tab.pair2.select[,4:6], 1, function(x) all(x == set.ind.C$names.statuses.tab[stAt2,])))
        initcond2.2[stAt2] <- sum(OUT2.P1[lastrow2.1, pos.to.calc2])
      }
      
     #calculate probabilities until C was sampled
      CoarseDeno <- 10
      OUT2.P2 <- ode(y = initcond2.2, times = seq(tAp,tCp, by = (tCp-tAp)/CoarseDeno), func = solveKFE, parms = list(model.mat=mat.ind.C), method = "radau")
      lastrow2.2 <- nrow(OUT2.P2)
      #check for vals < 0
      while(any(OUT2.P2[lastrow2.2,-1] < 0)){
        if(any(abs(OUT2.P2[lastrow2.2,which(OUT2.P2[lastrow2.2,-1] < 0)+1]) < 10^-30) | CoarseDeno > 2000){
          OUT2.P2[lastrow2.2,which(OUT2.P2[lastrow2.2,-1] < 0)+1] <- 0
        }else{
          CoarseDeno <- CoarseDeno*2
          OUT2.P2 <- ode(y = initcond2.2, times = seq(tAp,tCp, by = (tCp-tAp)/CoarseDeno), func = solveKFE, parms = list(model.mat=mat.ind.C), method = "radau")
          lastrow2.2 <- nrow(OUT2.P2)
          #print(CoarseDeno)
        }
      }
   
      #store.results
      p.c2.Tstart[tstart] <- sum(OUT2.P2[lastrow2.2,
                                         which(as.numeric(set.ind.C$names.statuses.tab[,2]) + 
                                               as.numeric(set.ind.C$names.statuses.tab[,3]) == div.C) + 1])
      
    #######################################
    }else if(tCp < tAp){#C SAMPLED BEFORE A
      #I already made sure that t0 != tAp
      if(t0 < tCp){
        #calculate probabilities until A was sampled at tAp
        CoarseDeno <- 10
        OUT2.P1T <- ode(y = initcond2.1, times = seq(t0,tCp,by=(tCp-t0)/CoarseDeno), func = solveKFE, parms = list(model.mat=mat.pair2), method = "radau")
        lastrow2.1 <- nrow(OUT2.P1T)
        #check for vals < 0
        while(any(OUT2.P1T[lastrow2.1,-1] < 0)){
          if(any(abs(OUT2.P1T[lastrow2.1,which(OUT2.P1T[lastrow2.1,-1] < 0)+1]) < 10^-20) | CoarseDeno > 2000){
            OUT2.P1T[lastrow2.1,which(OUT2.P1T[lastrow2.1,-1] < 0)+1] <- 0
          }else{
            CoarseDeno <- CoarseDeno*2
            OUT2.P1T <- ode(y = initcond2.1, times = seq(t0,tCp,by=(tCp-t0)/CoarseDeno), func = solveKFE, parms = list(model.mat=mat.pair2), method = "radau")
            lastrow2.1 <- nrow(OUT2.P1T)
            #print(CoarseDeno)
          }
        }
        ##########
        #chain 2 -> A infected C
        OUT2.P1 <- OUT2.P1T[,which(as.numeric(set.pair2$names.statuses.tab[,5]) + 
                                   as.numeric(set.pair2$names.statuses.tab[,6]) == div.C & 
                                   is.element(set.pair2$names.statuses.tab[,4], c("E","I"))) + 1]
      
        names.statuses.tab.pair2.select <- set.pair2$names.statuses.tab[which(as.numeric(set.pair2$names.statuses.tab[,5]) + 
                                                                                as.numeric(set.pair2$names.statuses.tab[,6]) == div.C & 
                                                                                is.element(set.pair2$names.statuses.tab[,4], c("E","I"))),]
        #set initial conditions for second step
        initcond2.2 <- numeric(set.ind.A$tot.statuses)
        for(stAt2 in 1:set.ind.A$tot.statuses){
          pos.to.calc2 <- which(apply(names.statuses.tab.pair2.select[,1:3], 1, function(x) all(x == set.ind.A$names.statuses.tab[stAt2,])))
          initcond2.2[stAt2] <- sum(OUT2.P1[lastrow2.1, pos.to.calc2])
        }
        #calculate probabilities until C was sampled
        CoarseDeno <- 10
        OUT2.P2 <- ode(y = initcond2.2, times = seq(tCp,tAp, by = (tAp-tCp)/CoarseDeno), func = solveKFE, parms = list(model.mat=mat.ind.A), method = "radau")
        lastrow2.2 <- nrow(OUT2.P2)
        #check for vals < 0
        while(any(OUT2.P2[lastrow2.2,-1] < 0)){
          if(any(abs(OUT2.P2[lastrow2.2,which(OUT2.P2[lastrow2.2,-1] < 0)+1]) < 10^-30) | CoarseDeno > 2000){
            OUT2.P2[lastrow2.2,which(OUT2.P2[lastrow2.2,-1] < 0)+1] <- 0
          }else{
            CoarseDeno <- CoarseDeno*2
            OUT2.P2 <- ode(y = initcond2.2, times = seq(tCp,tAp, by = (tAp-tCp)/CoarseDeno), func = solveKFE, parms = list(model.mat=mat.ind.A), method = "radau")
            lastrow2.2 <- nrow(OUT2.P2)            
          }
        }
      
        #store.results
        p.c2.Tstart[tstart] <- sum(OUT2.P2[lastrow2.2,
                                           which(as.numeric(set.ind.A$names.statuses.tab[,2]) == 0 &  
                                                 as.numeric(set.ind.A$names.statuses.tab[,3]) == div.A) + 1])
         }else if(t0 >= tCp){
            p.c2.Tstart[tstart] <- 0
         }
    #######################################
    }else if(tCp == tAp){#C SAMPLED BEFORE A
      
        #calculate probabilities until A was sampled at tAp
        CoarseDeno <- 10
        OUT2.P1T <- ode(y = initcond2.1, times = seq(t0,tCp,by=(tCp-t0)/CoarseDeno), func = solveKFE, parms = list(model.mat=mat.pair2), method = "radau")
        lastrow2.1 <- nrow(OUT2.P1T)
                #check for vals < 0
        while(any(OUT2.P1T[lastrow2.1,-1] < 0)){
          if(any(abs(OUT2.P1T[lastrow2.1,which(OUT2.P1T[lastrow2.1,-1] < 0)+1]) < 10^-30) | CoarseDeno > 2000){
            OUT2.P1T[lastrow2.1,which(OUT2.P1T[lastrow2.1,-1] < 0)+1] <- 0
          }else{
            CoarseDeno <- CoarseDeno*2
            OUT2.P1T <- ode(y = initcond2.1, times = seq(t0,tCp,by=(tCp-t0)/CoarseDeno), func = solveKFE, parms = list(model.mat=mat.pair2), method = "radau")
            lastrow2.1 <- nrow(OUT2.P1T)
            #print(CoarseDeno)
          }
        }
#        rowSums(OUT2.P1T[,-1])
        ##########
        #chain 2 -> A infected C
        #no second step as they are sampled at the same time
        #store.results
        p.c2.Tstart[tstart] <- sum(OUT2.P1T[lastrow2.1,
                                            which(set.pair2$names.statuses.tab[,3] == div.A &
                                                  set.pair2$names.statuses.tab[,6] == div.C &
                                                  set.pair2$names.statuses.tab[,1] == "I") + 1])
        
      }
      #cat(paste("\r", format(100*tstart/t0.len,nsmall = 1), "% done"))
    }
    p.c2.Tstart.allPairs[,Num] <- p.c2.Tstart
  
  #print(Num)
}

return(list(chain.AC.prob = p.c2.Tstart.allPairs, t0 = t0.mat))
}

##############################################
