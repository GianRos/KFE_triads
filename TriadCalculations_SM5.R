###################################################################
#Code to set up triads calculations (Section S5)
#Include codes to generate KFE for SEI model
#
#Date:
#07/05/2019 - First script
#10/05/2019 - Update
#25/05/2019 - Edits
#
#Author:
#Gianluigi Rossi
###################################################################
rm(list = ls())

#load code for pair calculations
source("190424_KFE_pairs_calc_noPrints.R")
#load code for triad calculations
source("190510_KFE_triad_calc_noPrints.R")

#check folder
getwd()
###################################################################
#prepare parameters
betaTB.est <- c(5.236895e-05, 7.711509e-05, 1.062255e-04)
latperTB.est <- c(12.88906,  472.52607, 1737.05323)
subInfTB.est <- c(5.813327e-06, 4.308351e-04, 3.924632e-03)

###################################################################
#create pairs.list ad hoc for triad calculation
life.exp <- 365*8
max.time <- 365*10
deltaT.vect <- seq(0, max.time, by = 365)
deltaT.vect.l <- length(deltaT.vect)

#Note: setting
maxSnps <- 1

pairs.list.fortri <- NULL

if(maxSnps > 0){
  for(S1 in 0:maxSnps){
    for(S2 in 0:maxSnps){
      pairs.list.fortri <- rbind(pairs.list.fortri,
                               data.frame(
                                 ind1 = "A",
                                 ind2 = "C",
                                 sample.t1 = 0,
                                 sample.t2 = deltaT.vect,
                                 snps.d1 =   S1,
                                 snps.d2 =   S2,
                                 stringsAsFactors = F
                               )
                              )
    }
  }
}else{
  pairs.list.fortri <- data.frame(
    ind1 = "A",
    ind2 = "C",
    sample.t1 = 0,
    sample.t2 = deltaT.vect,
    snps.d1 =   0,
    snps.d2 =   0,
    stringsAsFactors = F
  )
}


###################################################################

pairs.list.fortri.final <- NULL

for(SNPS in 0:maxSnps){
  if(SNPS == 0){
  pairs.list.fortri.cut <- subset(pairs.list.fortri, 
                                snps.d1 <= 0 & snps.d2 <= 0)
  }else{
    pairs.list.fortri.cut <- subset(pairs.list.fortri, 
                                    (snps.d1 > (SNPS-1) | snps.d2 > (SNPS-1)) & snps.d1 <= SNPS & snps.d2 <= SNPS)
  }

    part.calc.ABC <- KFE_SEI_triads_calc(pairs.list = pairs.list.fortri.cut,
                    all.Params = c("tau" = betaTB.est[2],
                                   "sigma" = 1/latperTB.est[2],
                                   "max.life.exp" = life.exp,
                                   "sub.inf" = subInfTB.est[2],
                                   "sub.lat" = 1), #
                    #starting time
                    t0.fixed = NA,
                    #number of initial time steps
                    integral.time = 100)
  print("ABC_done")
  part.calc.AC <- KFE_SEI_pairs_calc(pairs.list = pairs.list.fortri.cut,
                                     all.Params = c("tau" = betaTB.est[2],
                                                    "sigma" = 1/latperTB.est[2],
                                                    "max.life.exp" = life.exp,
                                                    "sub.inf" = subInfTB.est[2],
                                                    "sub.lat" = 1), #
                                     #starting time
                                     t0.fixed = NA,
                                     #number of initial time steps
                                     integral.time = 100)
  print("AC_done")
  
  pairs.list.fortri.cut$ProbABC <- apply(part.calc.ABC$chain.AC.prob, 2, max)
  pairs.list.fortri.cut$ProbAC <- apply(part.calc.AC$chain.AC.prob, 2, max)

  pairs.list.fortri.final <- rbind(pairs.list.fortri.final, pairs.list.fortri.cut)
  
  save.image("KFE_SEI_triads_calculations.Rdata")
  print(SNPS)
}



###################################################################




