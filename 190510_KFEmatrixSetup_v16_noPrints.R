#rm(list = ls())
############################################
#Scope: Function to create a matrix for the Folmogorov Forward Equation model
#Date first finalised: 25/01/2019

#Notes/comments: 
#1.this code was built for a SETI bovine TB model. For different type of models edits are needed. 
#2.we assume that the infection start with the first individualand only he can infect other ones. To permit other transmission routes modifications are needed.
#3.tested against ODE solver with code 180128_KFEmatrixsetup_and_check.R, 2 individuals 1 SNP subs SETI model

#Version:
#v000-180126 -> base
#v01-180209 -> one input more: which chain (ex. individual 1 to 2, 2 to 3, 1 to 3, etc)
#v02-180214 -> no matrix expoenntial in output, too slow for big matrices
#v03-180216 -> the second SNPs index for newly infected start always from 0
#           -> the index individual does not show SNPs substitutions before infecting 
#v04-180222 -> divide former function argument SIMPLIFY.MAT in three argments:
#              1. index.snps, if FALSE I do not record the SNPs substitution in index individual before infecting,only status I00 can infect
#              2. index.ind, if FALSE the first individual is not the index individual, thus I do not remove the S status, as well as all status with the second SNPs > 0
#              3. simplify.mat
#v05-180410 -> changed inputs
#              added chain.start, TRUE if the chain I am analysing starts here, FALSE otherwise
#              added second.inf, TRUE if the first individual in the matrix was the second in a chain (remove useless status)
#              removed index.ind
#              added rem.idx.snps, TRUE if I want to reliminate the possibility of SNPs before first infection from index ind (#1) to the 2nd
#              removed index.snps
#              removed transm.pair : now it works as a chain, 1 infect 2, which infect, 3, etc.
#              added life.exp: life expectancy parameter (rate death.rate)
# NB: all SNPs of an individual are tracked from the previous one in the chain. If I want to track the SNPs from 2 or more chain jumps (ex indi 3 w/respect to A) I have to compute manually after
#v06-180417 -> changes in how to consider the absorbing state:
#              from now on, all SNPs appearances over sub.Max will follor in the R category (removed), not tracked, but they will still happen#           -> changed inputs
#              added ind.to.track: vector of number who identfy the individuals for which I want to exactly track theri SNPs, so after they reached the maximum number allowed (sub.max) of SNPs the are REMOVED
#v07-180420 -> Added explicit R status
#v08-180426 -> no Test positive class
#v09-180510 -> no removal for Exposed
#v10-180516 -> different removal rate for Exposed and Infectious
#v11-190125 -> removed status R to simplify matrix
#v11_noPrints-190311 -> removed message printing
#v12_noPrints-190405 -> check installation iterpc package outside function
#v13_noPrints-190413 -> restyled matrix and status table definition
#                    -> addition of a D status for deaths
#                    -> BIG CHANGE: ind.to.track replaced with unsampled.individuals (FUNCTION INPUT, TO FIX IN FOLLOWING FILES)
#                    -> BIG CHANGE: the function becomes specific for each type of model (at the moment SEI, will follow SIR, SI, SEIR), no need to specify model
#                    -> iterpc package NOT USED anymore
#                    -> When individuals reach maxmimum SNPs permitted they stay in same class
#v14_noPrints-190420 -> removed all.alive, its functions taken over by unspampled.individuals
#v15_noPrints-190420 -> if n.ind > 1 initial susceptibles don't die
#                    -> life.exp substituteb with max.life.exp
#v16_noPrints-190510 -> 

#Author:
#Gianluigi Rossi
############################################
#Packages/libraries required:
require(deSolve)

############################################
#Input: 
#model type, #individuals, #maximum SNPs substitutions allowed, model parameters, reduce the number of status to the necessary (yes or no)

#Outputs: 
#matrix, exponential matrix, number of status, names of status

############################################
#start function to create matrix specific for SEI models 
set.KFE.system.SEI <- function(n.ind, #number of individuals in the system
                           sub.max, #maximum of SNPs  allowed in the system
                           params, #vector with the parameters (tau, sigma, gamma, death.rate, life.exp, sub.inf, sub.lat)
                           unsampled.individuals = NA, #individuals in the chain but not sampled, so I do not track the number of SNPs in the results
                           simplify.mat = TRUE, #logical, if true I reduce the matrix size by reducing ther number of status
                           #the following arguments works only if simplify.mat = TRUE
                           div.snps.only = TRUE, #no SNPs substitutions before infection
                           chain.start = TRUE, #is the first step of a chain?
                           second.inf = FALSE) #is the first individual the second of a chain?
    with(as.list(params),{
    
    all.compart <- c("S","E","I")
      
    ###temporary
    #n.ind <- 3
    #sub.max <- 1
    ###end
    #print(params)
    
    death.rate <- -log(0.001)/max.life.exp
    
    n.compart <- length(all.compart)  
    ind.tot.statuses <- 2 + (n.compart-1)*(sub.max+1)*(sub.max+1)
     
    ind.SNPs.comb <- expand.grid(0:sub.max,0:sub.max)

    ind.statuses.table <- matrix(-1, ncol = 3, nrow = ind.tot.statuses)
    rownames(ind.statuses.table) <- paste("Status", 1:nrow(ind.statuses.table), sep ="")
    colnames(ind.statuses.table) <- c("Disease","NumHeredSNPs","NumDivSNPs")

    #single individual starts with susceptible status
    ind.statuses.table[1,1] <- "S"
    #status E
    ind.statuses.table[1:((sub.max+1)*(sub.max+1)) + 1,1] <- "E"
    ind.statuses.table[1:((sub.max+1)*(sub.max+1)) + 1,2] <- ind.SNPs.comb$Var1
    ind.statuses.table[1:((sub.max+1)*(sub.max+1)) + 1,3] <- ind.SNPs.comb$Var2
    #status I
    ind.statuses.table[1:((sub.max+1)*(sub.max+1)) + ((sub.max+1)*(sub.max+1)) + 1,1] <- "I"
    ind.statuses.table[1:((sub.max+1)*(sub.max+1)) + ((sub.max+1)*(sub.max+1)) + 1,2] <- ind.SNPs.comb$Var1
    ind.statuses.table[1:((sub.max+1)*(sub.max+1)) + ((sub.max+1)*(sub.max+1)) + 1,3] <- ind.SNPs.comb$Var2
    #last status is death
    ind.statuses.table[nrow(ind.statuses.table),1] <- "D"
  
    #Number of all potential statuses
    all.statuses <- ind.tot.statuses^n.ind
    
    #create all potential statuses table
    if(n.ind > 1){
      step1str <- paste(c("1:ind.tot.statuses",rep(",1:ind.tot.statuses",(n.ind-1))), collapse = "")
      posi.vect <- eval(parse(text = paste("expand.grid(",step1str,")",sep = "")))
      step2str <- paste(c("cbind(ind.statuses.table[posi.vect[,1],]",rep(paste(",ind.statuses.table[posi.vect[,",2:n.ind,"],]")), ")"),collapse = "")
      all.statuses.table <- eval(parse(text = step2str))
      rownames(all.statuses.table) <- paste("Status", 1:nrow(all.statuses.table))
      colnames(all.statuses.table) <- paste(rep(colnames(ind.statuses.table),n.ind), rep(1:n.ind,each=3),sep = "")
    }else{
      all.statuses.table <- ind.statuses.table
    }

    #Matrix names: combination of the 3 statuses
    #names.statuses <- apply(table.statuses, 1, paste, collapse = "")

    #Before building the matrix I remove the status that are not useful, following assumptions (listed below)
    status.to.remove <- vector()
    #Simplify only if chains are longer than 1 individual
    if(simplify.mat & n.ind > 1){

      #all.statuses.table[status.to.remove,]
      #status.to.remove <- vector()
      
        #If the chain starts here (chain.start = TRUE)
        #NB: the matrix can be computed for the second part of the chain, in which the computation picks up the result of the previous one
        if(chain.start){
            #1. I remove the index individual (#1) S status
            status.to.remove <- c(status.to.remove, which(all.statuses.table[,1] == "S"))
            #2. Index individual (1) has no hereditary SNPs
            if(sub.max>0){
                status.to.remove <- c(status.to.remove, which(as.numeric(all.statuses.table[,2]) > 0))
            }
        
            #If I do not want to track the SNPs appearing on index individual before infecting others (tracking only the divergent ones)
            #If the are more than zero SNPs allowed
            if(div.snps.only & sub.max > 0){
            #3. Index individual (1) SNPs substitutions tracked after infection, so removed NumHerSNPs > 0 for second individual
              status.to.remove <- c(status.to.remove, which(as.numeric(all.statuses.table[,5]) > 0))
              status.to.remove <- c(status.to.remove, which(as.numeric(all.statuses.table[,3]) > 0 & all.statuses.table[,4] == "S"))
            }

        }else if(second.inf & div.snps.only & sub.max > 0){
            status.to.remove <- c(status.to.remove, which(as.numeric(all.statuses.table[,2]) > 0))
        }
        
        #Individuals can only be infected by the previous one in the chain (1>2, 2>3, 3>4, etc.)
        for(ind in 2:n.ind){
          status.to.remove <- c(status.to.remove, which(!is.element(all.statuses.table[,((ind-1)*3 - 2)], c("I","D")) & all.statuses.table[,(ind*3 - 2)] != "S"))
        }
      
        if(is.na(unsampled.individuals)){
          #if not interested in statuses in which one individual is dead
          remove.dead <- which(apply(all.statuses.table, 1, function(x) any(x == "D")))
          #keep the status in which everyone is dead to use as a destination for deaths
          remove.dead <- remove.dead[-which(remove.dead == which(apply(all.statuses.table[,seq(1, (n.ind*3 - 2), by = 3)], 1, function(x) all(x == rep("D",n.ind)))))]
          status.to.remove <- c(status.to.remove, remove.dead)
        }else{
          #we keep the status in which the not sampled individuals are dead, because we do not care to find them alive
          if(n.ind == 2){
              remove.dead <- which(all.statuses.table[,(1:n.ind)[-unsampled.individuals] * 3 - 2] == "D")
          }else if(n.ind >2){
              remove.dead <- which(apply(all.statuses.table[,(1:n.ind)[-unsampled.individuals] * 3 - 2], 1, function(x) any(x == "D")))
          }
          #keep the status in which everyone is dead to use as a destination for deaths
          remove.dead <- remove.dead[-which(remove.dead == which(apply(all.statuses.table[,seq(1, (n.ind*3 - 2), by = 3)], 1, function(x) all(x == rep("D",n.ind)))))]
          status.to.remove <- c(status.to.remove, remove.dead)
        }
    }else if(simplify.mat){
      if(chain.start){
        #1. I remove the index individual (#1) S status
        status.to.remove <- c(status.to.remove, which(all.statuses.table[,1] == "S"))
        #2. Index individual (1) has no hereditary SNPs
        if(sub.max>0){
          status.to.remove <- c(status.to.remove, which(as.numeric(all.statuses.table[,2]) > 0))
        }
      }else if(second.inf){
        status.to.remove <- c(status.to.remove, which(as.numeric(all.statuses.table[,2]) > 0))
      }
    }
    #remove doubles
    status.to.remove <- unique(status.to.remove)
    #remove status
    if(length(status.to.remove)>0){
      all.statuses.table <- all.statuses.table[-status.to.remove,]
    }
    #update all.statuses.table characteristics
    tot.statuses <- nrow(all.statuses.table)
    row.names(all.statuses.table) <- paste("Status", 1:tot.statuses, sep = "")
    names.statuses <- row.names(all.statuses.table)
    
    #Matrix
    Q.mat <- matrix(0,
                    nrow = tot.statuses,# for S, 3 for S, E, T x number of substitution allowed
                    ncol = tot.statuses)
    colnames(Q.mat) <- names.statuses
    rownames(Q.mat) <- names.statuses
  
    #############

    for(i in 1:tot.statuses){
        #calculate rates per individuals i <- 1
        for(NI in 1:n.ind){
            #identify epidemiological status NI <- 1
            curr.epi.stat <- all.statuses.table[i,NI*3-2]       
        
            #select the potential destination statuses, based on status of other individuals (if.n.ind > 1)
            if(n.ind ==  1){
                j.select.1 <- 1:tot.statuses
            }else{
                j.select.0 <- which(apply(all.statuses.table[,-c(NI*3-(2:0))], 1, function(X) all(X == all.statuses.table[i,-c(NI*3-(2:0))])))
                j.select.1 <- j.select.0[-which(j.select.0 == i)]
            }
            
            #INFECTION -> S to E0 
            if(curr.epi.stat == "S" & n.ind > 1 & NI > 1){

                #NB: transmissions always happens in scale (1>2, 2>3, 3>4, etc...NOT 1>3)
                transm.ind <- NI-1
                if(all.statuses.table[i,transm.ind*3-2] == "I"){
                    curr.SNP.infect <- all.statuses.table[i,transm.ind*3]
                                        #infection are passed with substituted SNPs
                    j.select.2 <- j.select.1[which(all.statuses.table[j.select.1,NI*3-2] == "E" &
                                                   all.statuses.table[j.select.1,NI*3-1] == curr.SNP.infect & 
                                                   all.statuses.table[j.select.1,NI*3] == 0)]
                    Q.mat[i,j.select.2] <- Q.mat[i,j.select.2]+tau
                }
                
                #DEATH -> S to D
                #j.select.2 <- which(all.statuses.table[,NI*3-2] == "D")
                #Q.mat[i,j.select.2] <- Q.mat[i,j.select.2] + death.rate
     
            }else{
                orig.SNP.ind <- all.statuses.table[i,NI*3-1]
                curr.SNP.ind <- all.statuses.table[i,NI*3]
                
                #END LATENT PERIOD -> EXX to IXX
                if(curr.epi.stat == "E"){
                    j.select.2 <- j.select.1[which(all.statuses.table[j.select.1,NI*3-2] == "I" &
                                                   all.statuses.table[j.select.1,NI*3-1] == orig.SNP.ind &
                                                   all.statuses.table[j.select.1,NI*3] == curr.SNP.ind)]
                    Q.mat[i,j.select.2] <- sigma
                }
          
              
                #SNPs SUBSTITUTION
                if(curr.SNP.ind < sub.max){
                    j.select.2 <- j.select.1[which(all.statuses.table[j.select.1,NI*3-2] == curr.epi.stat &
                                                     all.statuses.table[j.select.1,NI*3-1] == orig.SNP.ind &
                                                     all.statuses.table[j.select.1,NI*3] == (as.numeric(curr.SNP.ind)+1))]    
                    if(curr.epi.stat == "I"){
                        Q.mat[i,j.select.2] <- sub.inf
                    }else{
                        Q.mat[i,j.select.2] <- sub.lat*sub.inf
                    }
                }
                #individuals with more SNPs than maximum go to Death
                #NB: this might be changed later on
                 #else if(curr.SNP.ind == sub.max){
                    #j.select.2 <- j.select.1[which(substr(table.statuses[j.select.1,NI],1,1) == "R")]
                    #if(curr.epi.stat == "I"){
                    #    Q.mat[i,j.select.2] <- sub.inf
                    #}else{
                    #    Q.mat[i,j.select.2] <- sub.lat*sub.inf
                    #}
                #}

                #REMOVAL -> IXX to R
                #if(curr.epi.stat == "I"){
                #    j.select.2 <- j.select.1[which(all.statuses.table[j.select.1,NI*3-2] == "R")]
                #    Q.mat[i,j.select.2] <- Q.mat[i,j.select.2] + death.rate
                #}
                
                #DEATH -> any to D
                if(curr.epi.stat != "D"){
                  if(!is.element(NI, unsampled.individuals)){
                    #j.select.2 <- j.select.1[which(all.statuses.table[j.select.1,NI*3-2] == "D")]
                    #Q.mat[i,j.select.2] <- Q.mat[i,j.select.2] + death.rate
                    #}else{
                    j.select.2 <- which(all.statuses.table[,NI*3-2] == "D")
                    Q.mat[i,j.select.2] <- Q.mat[i,j.select.2] + death.rate      
                  }
                }
              }
        }
        Q.mat[i,i] <- -sum(Q.mat[i,-i])
       # cat(paste("\r status", i, "of", tot.statuses))
    }
    #all.statuses.table[j.select.2,]

    COUNT <- 0
    while(sum(Q.mat) < 0 & COUNT < 50){
      Q.mat[which(rowSums(Q.mat)>0),ncol(Q.mat)] <- Q.mat[which(rowSums(Q.mat)>0),ncol(Q.mat)] - rowSums(Q.mat)[which(rowSums(Q.mat)>0)]
      COUNT <- COUNT +1
      #print(COUNT)
    }
    if(any(Q.mat[,ncol(Q.mat)] < 0)){
      Q.mat[which(Q.mat[,ncol(Q.mat)] < 0),ncol(Q.mat)] <- 0
    } 
    
    return(list(
        KFE.matrix=Q.mat, #Kolmogorov Forward Equations system matrix
        names.statuses.tab = all.statuses.table, #table with statuses, each column correspond to an individual
        tot.statuses = tot.statuses #final number of statuses
    ))
    })
############################################
