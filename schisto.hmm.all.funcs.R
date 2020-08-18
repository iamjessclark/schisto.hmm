# Functions for the models and data manipulation for 
# A hidden markov model clearance and reinfection, schisto
# Dr J. Clark, Dr P. H. L. Lamberton, Dr J. Prada

#### MODEL FUNCTIONS ####
#### do.kk.run ####

# this is just the function for the KK data

do.KK.run <- function (N, Ti, R, KK) {
  ## Set seed ##
  .RNG.seed <- function(chain)
    return( switch(chain, "1"= 1, "2"= 2) )
  .RNG.name <- function(chain)
    return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )
  
  ## Initialize Status ##
  IniStatus <- rep(1,N)
  Clearance <- matrix(0,nrow=N,ncol=Ti-1)
  Reinfec <- matrix(1,nrow=N,ncol=Ti-1)
  sigma<- c(rep(NA,3),0.01)
  
  m <- "model {
    
    # Prior prevalence / clearance / reinfection #

    prev ~ dbeta(1,1)
    
    clearance[2] ~ dbeta(1,1)
    clearance[3] ~ dbeta(1,1)
    clearance[4] ~ dbeta(1,1)
    
    reinfec[2] <- 0
    reinfec[3] ~ dbeta(1,1)
    reinfec[4] ~ dbeta(1,1)
    
    # Prior KKs #
    
    rtnb ~ dgamma(286.3235, 223.4650)
    sh ~ dgamma(83.88592,125.70331)
    rt1 ~ dbeta(47.13542,633.08366)
    rt <- rt1/(1-rt1)
    
    # accounting for the difference in length of time steps because there will be different amounts of variance
    # these are the priors for the random walk 
    sigma[4] ~ dgamma(0.001, 0.001) # sigma when the transition from timesteps up to 4
    sigma[3] <- sigma[4]*6/15 # sigma when the transition from timestep 2 up to 3 scaling down the variance from the previous step
    sigma[2] <- sigma[4]*3/15 # sigma when the transition from timestep 
    sigma[1] <- 0 #not used but needed to initialize sigma
    
    # model 
  for(n in 1:N){	# run through pop
    
    IniStatus[n] ~ dbern(prev) # the probability of you being infected at the beginning is dependent on the prevalence of infection, doesn't this account for some of them being infected and some of them not at baseline because some zero counts might still be infected at low level infections and so might show up as a 1 or 2 on the CCA scale 
    Status[n,1] <- IniStatus[n] # individual n's status at t1 is taken from the draw above
    
    # this are the true baseline counts set outside of the time loop for the auto-regressive part 
    tKK[n,1,1] <- 0
    tKK[n,2,1] ~ dgamma(sh,rt) 
    
    for (t in 1:Ti){ # run through timepoints from timept 2 onwards
      lambda[n,t] <- tKK[n,Status[n,t]+1,t] # this needs to go here because each time the r loop reiterates it will replace this value
        for(r in 1:R){ # run through repeated measures to set the baseline KK over the repeated measures, 
           KK[n,t,r] ~ dnegbin(rtnb/(lambda[n,t]+rtnb),rtnb)
        } # end or r loop
     } # end timestep 1 t loop
    
    for (t in 2:Ti){ # HMM component
      tKK[n,1,t] <- 0
      tKK[n,2,t] ~ dnorm(tKK[n,Status[n,t-1]+1,t-1],sigma[t]) # taking the mean in the previous time step + some sort of movement around that mean, the reason it is a normal dist is because this is the movement around the previous mean
    
      Clearance[n,t-1] ~ dbern(clearance[t])
      Reinfec[n,t-1] ~ dbern(reinfec[t])
      Status[n,t] <- Status[n,t-1] * (1-Clearance[n,t-1]) + (1-Status[n,t-1]) * Reinfec[n,t-1]
      
      switch[n,t] <- abs(Status[n,t]-Status[n,t-1])      # this is a dummy variable looking at how many times an individual switches status over the run

      } # end hmm time
  } # end n loop
    
    #inits# .RNG.seed, .RNG.name, IniStatus, Clearance, Reinfec   
    #data# N, Ti, R, KK
    #monitor# prev, clearance, reinfec, rtnb, sh, rt, Status
}"
  
  
  # Run model #
  Results <- run.jags(m, burnin=10000, sample=20000, n.chains=2, jags.refresh = 1, method = 'parallel',
                      plots = F, silent.jags = F)
  return(Results)
}

#### do.run ####

# This second function has both CCA and KK in it 

do.run <- function (N, Ti, R, KK, CCA) {
  ## Set seed ##
  ##this is so that they start on different values ##
  .RNG.seed <- function(chain)
    return( switch(chain, "1"= 1, "2"= 2) )
  .RNG.name <- function(chain)
    return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )
  
  ## Initialize Status ##
  IniStatus <- rep(1,N)
  Clearance <- matrix(0,nrow=N,ncol=Ti-1)
  Reinfec <- matrix(1,nrow=N,ncol=Ti-1)
  prob <- array(rep(CCA,2),dim=c(N,Ti,2))
  sigma<- c(rep(NA,3),0.01)
  
  m <- "model {
  # Prior prevalence / clearance / reinfection #
  prev ~ dbeta(1,1)

  clearance[2] ~ dbeta(1,1)
  clearance[3] ~ dbeta(1,1)
  clearance[4] ~ dbeta(1,1)
  
  reinfec[2] <- 0
  reinfec[3] ~ dbeta(1,1)
  reinfec[4] ~ dbeta(1,1)
  
  # Prior KKs 
  rtnb ~ dgamma(286.3235, 223.4650)
  sh ~ dgamma(83.88592,125.70331)
  rt1 ~ dbeta(47.13542,633.08366)
  rt <- rt1/(1-rt1)
  
  # Priors for the random walk #
  # accounting for the difference in length of time steps because there will be different amounts of variance
  sigma[4] ~ dgamma(0.001, 0.001) # sigma when the transition from timesteps up to 4
  sigma[3] <- sigma[4]*6/15 # sigma when the transition from timestep 2 up to 3 scaling down the variance from the previous step
  sigma[2] <- sigma[4]*3/15 # sigma when the transition from timestep 
  sigma[1] <- 0 #not used but needed to initialize sigma
  
  # Priors CCA #
  #k ~ dnorm(0.063291330, 1/0.009817479^2)
  #intercept ~ dunif(0.0139464, 8.5045100)
  
  # I have changed these from the informative priors because the logistic curve is now not the same as before
  k ~ dgamma(0.001, 0.001)
  intercept ~ dgamma(0.001, 0.001)
  
  # MODEL 
  for(n in 1:N){	# run through pop

    IniStatus[n] ~ dbern(prev)
    Status[n,1] <- IniStatus[n]
    
    # these are the true baseline counts set outside of the time loop for the auto-regressive part 
    tKK[n,1,1] <- 0
    tKK[n,2,1] ~ dgamma(sh,rt) 
 
    prob[n,1,1] ~ dnorm(0,3.093451)
    prob[n,1,2] ~ dnorm(4 / (1 + exp(-k*(tKK[n,2,1]-intercept))),3.093451)

    CCA[n,1] ~ dround(prob[n,1,Status[n,1]+1],0)
    
    for (t in 1:Ti){ # run through timepoints
      lambda[n,t] <- tKK[n,Status[n,t]+1,t] # this needs to go here because each time the r loop reiterates it will replace this value
        for( r in 1:R){  # run through repeat measures
          KK[n,t,r] ~ dnegbin(rtnb/(lambda[n,t]+rtnb),rtnb) # generating the data with noise and then sampling from the dataset with the gamma sh/rt parameters?
      } # end or r loop
    } # end timestep 1 t loop
    
    for (t in 2:Ti){ # HMM component
      tKK[n,1,t] <- 0                                                                                                                                                                                                                                                                                               
      tKK[n,2,t] ~ dnorm(tKK[n,Status[n,t-1]+1,t-1],sigma[t])
      
      prob[n,t,1] ~ dnorm(0,3.093451)
      prob[n,t,2] ~ dnorm(4 / (1 + exp(-k*(tKK[n,2,t]-intercept))),3.093451)
       
      CCA[n,t] ~ dround(prob[n,t,Status[n,t]+1],0)
      
      Clearance[n,t-1] ~ dbern(clearance[t])
      Reinfec[n,t-1] ~ dbern(reinfec[t])

      Status[n,t] <- Status[n,t-1] * (1-Clearance[n,t-1]) + (1-Status[n,t-1]) * Reinfec[n,t-1] 
    }
    
  }  
  
  #inits# .RNG.seed, .RNG.name, IniStatus, Clearance, Reinfec, sigma, prob   
  #data# N, Ti, R, KK, CCA
  #monitor#  prev, clearance, reinfec, rtnb, k, intercept, sh, rt, Status, CCA
}"
  
  # Run model #
  Results <- run.jags(m, burnin=10000, sample=20000, n.chains=2, jags.refresh = 1, method = 'parallel',
                      plots = F, silent.jags = F)
  return(Results)
}

#### a function for the 10 step scale on CCA diagnostics ####
do.10.run <- function (N, Ti, R, KK, CCA10) {
  ## Set seed ##
  ##this is so that they start on different values ##
  .RNG.seed <- function(chain)
    return( switch(chain, "1"= 1, "2"= 2) )
  .RNG.name <- function(chain)
    return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )
  
  ## Initialize Status ##
  IniStatus <- rep(1,N)
  Clearance <- matrix(0,nrow=N,ncol=Ti-1)
  Reinfec <- matrix(1,nrow=N,ncol=Ti-1)
  prob <- array(rep(CCA10,2),dim=c(N,Ti,2))
  sigma<- c(rep(NA,3),0.01)
  
  m <- "model {
  
  # MODEL 
   # Prior prevalence / clearance / reinfection #
  prev ~ dbeta(1,1)
  
  clearance[2] ~ dbeta(1,1)
  clearance[3] ~ dbeta(1,1)
  clearance[4] ~ dbeta(1,1)
  
  # reinfection at the first time step set to 0 to reduce the number of switches
  reinfec[2] <- 0
  reinfec[3] ~ dbeta(1,1)
  reinfec[4] ~ dbeta(1,1)
  
  # Prior KKs 
  rtnb ~ dgamma(286.3235, 223.4650)
  sh ~ dgamma(83.88592,125.70331)
  rt1 ~ dbeta(47.13542,633.08366)
  rt <- rt1/(1-rt1)
  
  # Priors for the random walk #
  # accounting for the difference in length of time steps because there will be different amounts of variance
  sigma[4] ~ dgamma(0.001, 0.001) # sigma when the transition from timesteps up to 4
  sigma[3] <- sigma[4]*6/15 # sigma when the transition from timestep 2 up to 3 scaling down the variance from the previous step
  sigma[2] <- sigma[4]*3/15 # sigma when the transition from timestep 
  sigma[1] <- 0 #not used but needed to initialize sigma
  
  # Priors CCA gscore - not the same as the normal CCA so not using priors from the previous model #
  k ~ dgamma(0.001, 0.001)
  intercept ~ dgamma(0.001, 0.001)
  
  for(n in 1:N){	# run through pop
  
    IniStatus[n] ~ dbern(prev)
    Status[n,1] <- IniStatus[n]
    
    # these are the true baseline counts set outside of the time loop for the auto-regressive part 
    tKK[n,1,1] <- 0
    tKK[n,2,1] ~ dgamma(sh,rt) 

    prob[n,1,1] ~ dnorm(0,1.093606)
    prob[n,1,2] ~ dnorm(9 / (1 + exp(-k*(tKK[n,2,1]-intercept))),1.093606) 

    
    CCA10[n,1] ~ dround(prob[n,1,Status[n,1]+1],0)
    
    for (t in 1:Ti){ # run through timepoints
          lambda[n,t] <- tKK[n,Status[n,t]+1,t] 
          for( r in 1:R){  # run through repeat measures
            KK[n,t,r] ~ dnegbin(rtnb/(lambda[n,t]+rtnb),rtnb) 
            } # end or r loop
        } # end timestep 1 t loop
    
    for (t in 2:Ti){ # HMM component
          tKK[n,1,t] <- 0
          tKK[n,2,t] ~ dnorm(tKK[n,Status[n,t-1]+1,t-1],sigma[t]) 

          prob[n,t,1] ~ dnorm(0,1.093606) 
          prob[n,t,2] ~ dnorm(9 / (1 + exp(-k*(tKK[n,2,t]-intercept))),1.093606)
          
          CCA10[n,t] ~ dround(prob[n,t,Status[n,t]+1],0)
          
          Clearance[n,t-1] ~ dbern(clearance[t])
          Reinfec[n,t-1] ~ dbern(reinfec[t])
          
          Status[n,t] <- Status[n,t-1] * (1-Clearance[n,t-1]) + (1-Status[n,t-1]) * Reinfec[n,t-1] 
          
          switch[n,t] <- abs(Status[n,t]-Status[n,t-1])      # this is a dummy variable looking at how many times an individual switches status over the run
        } # t loop
  
  } # N loop   
  
 
  #inits# .RNG.seed, .RNG.name, IniStatus, Clearance, Reinfec, prob,  sigma
  #data# N, Ti, R, KK, CCA10
  #monitor#  prev, clearance, reinfec, rtnb, k, intercept, sh, rt, Status, prob   
  
}"
  
  # Run model #
  Results <- run.jags(m, burnin=10000, sample=20000, n.chains=2, jags.refresh = 1, method = 'parallel',
                      plots = F, silent.jags = F)
  return(Results)
}

# ---------------------------------------------------------------------------------------------------------------

# keeping for now becuase it took me so long to get the blasted thing to run ill be damned if I am deleting it. 
# but doesn't actually work, wildly underestimates prevalence 

do.10.run.lomax <- function (N, Ti, R, KK, CCA10) {
  ## Set seed ##
  ##this is so that they start on different values ##
  .RNG.seed <- function(chain)
    return( switch(chain, "1"= 1, "2"= 2) )
  .RNG.name <- function(chain)
    return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )
  
  ## Initialize Status ##
  IniStatus <- rep(1,N)
  Clearance <- matrix(0,nrow=N,ncol=Ti-1)
  Reinfec <- matrix(1,nrow=N,ncol=Ti-1)
  prob <- array(rep(CCA10,2),dim=c(N,Ti,2))
  tau<-10
  sigma<- c(rep(NA,3),0.01)

  m <- "model {
 # model here is written backwards because that is how rjags reads them (bottom to top) 
  # MODEL 
  for(n in 1:N){	# run through pop
  
    IniStatus[n] ~ dbern(prev)
    Status[n,1] <- IniStatus[n]
    
    # these are the true baseline counts set outside of the time loop for the auto-regressive part 
    tKK[n,1,1] <- 0
    tKK[n,2,1] ~ dgamma(sh,rt) 
  
    #prob[n,1,1] ~ dnorm(munorm, tau2)
    #prob[n,1,1] ~ dpois(muPo)T(,9)
    #prob[n,1,1] ~ dgamma(sh2, rate)
    #prob[n,1,1] ~ dlomax(alphaPAR, sigmaPAR)
    prob[n,1,1] ~ dpar1(alphaPAR, sigmaPAR)
    prob[n,1,2] ~ dnorm(9 / (1 + exp(-k*(tKK[n,2,1]-intercept))),tau) 
    
    CCA10[n,1] ~ dround(prob[n,1,Status[n,1]+1],0)
    
    for (t in 1:Ti){ # run through timepoints
          lambda[n,t] <- tKK[n,Status[n,t]+1,t] 
          for( r in 1:R){  # run through repeat measures
            KK[n,t,r] ~ dnegbin(rtnb/(lambda[n,t]+rtnb),rtnb) 
            } # end or r loop
        } # end timestep 1 t loop
    
    for (t in 2:Ti){ # HMM component
          tKK[n,1,t] <- 0
          tKK[n,2,t] ~ dnorm(tKK[n,Status[n,t-1]+1,t-1],sigma[t]) 
          
          #prob[n,t,1] ~ dnorm(munorm, tau2)
          #prob[n,t,1] ~ dpois(muPo)T(,9)
          #prob[n,t,1] ~ dgamma(sh2, rate)
          prob[n,t,1] ~ dlomax(alphaPAR, sigmaPAR) # support begins at 0
          #prob[n,t,1] ~ dpar1(alphaPAR, sigmaPAR) # doesn't support 0 so won't work
          prob[n,t,2] ~ dnorm(9 / (1 + exp(-k*(tKK[n,2,t]-intercept))),tau)
          
          CCA10[n,t] ~ dround(prob[n,t,Status[n,t]+1],0)
          
          Clearance[n,t-1] ~ dbern(clearance[t])
          Reinfec[n,t-1] ~ dbern(reinfec[t])
          
          Status[n,t] <- Status[n,t-1] * (1-Clearance[n,t-1]) + (1-Status[n,t-1]) * Reinfec[n,t-1] 
          
          switch[n,t] <- abs(Status[n,t]-Status[n,t-1])      # this is a dummy variable looking at how many times an individual switches status over the run
        } # t loop
  
  } # N loop   
  
  # Prior prevalence / clearance / reinfection #
  prev ~ dbeta(1,1)
  
  clearance[2] ~ dbeta(1,1)
  clearance[3] ~ dbeta(1,1)
  clearance[4] ~ dbeta(1,1)
  
  # reinfection at the first time step set to 0 to reduce the number of switches
  reinfec[2] <- 0
  reinfec[3] ~ dbeta(1,1)
  reinfec[4] ~ dbeta(1,1)
  
  # Prior KKs 
  rtnb ~ dgamma(286.3235, 223.4650)
  sh ~ dgamma(83.88592,125.70331)
  rt1 ~ dbeta(47.13542,633.08366)
  rt <- rt1/(1-rt1)
  
  # Priors for the random walk #
  # accounting for the difference in length of time steps because there will be different amounts of variance
  sigma[4] ~ dgamma(0.001, 0.001) # sigma when the transition from timesteps up to 4
  sigma[3] <- sigma[4]*6/15 # sigma when the transition from timestep 2 up to 3 scaling down the variance from the previous step
  sigma[2] <- sigma[4]*3/15 # sigma when the transition from timestep 
  sigma[1] <- 0 #not used but needed to initialize sigma
  
  # Priors CCA gscore - not the same as the normal CCA so not using priors from the previous model #
  k ~ dnorm(0, 1/10^10)
  intercept ~ dunif(0, 100)
  tau ~ dgamma(0.001,0.001)
  alphaPAR ~ dgamma(0.001, 0.001)
  sigmaPAR ~ dunif(0.000001, 10000)
  
  #inits# .RNG.seed, .RNG.name, IniStatus, Clearance, Reinfec, prob, tau, sigma
  #data# N, Ti, R, KK, CCA10
  #monitor#  prev, clearance, reinfec, rtnb, tau, k, intercept, sh, rt, Status
  
}"
  
  # Run model #
  Results <- run.jags(m, burnin=5000, sample=10000, n.chains=2, jags.refresh = 1, 
                      method = 'rjags', modules="runjags", plots = F, silent.jags = F)
  load.runjagsmodule()
  return(Results)
}
 #### DATA FUNCTIONS ####

#### Data loading ####

#### Kato Katz Data ####

loadKKdata <- function(nameFile){
  dt <- read.csv(nameFile) # read data file
  
  # name different weeks
  dt$dateN <- NA
  dt$dateN[which(dt$date=="25/09/2017" | dt$date=="26/09/2017" | dt$date=="27/09/2017"
                 | dt$date=="28/09/2017" | dt$date=="29/09/2017" | dt$date=="02/10/2017")] <- "Pre-T"
  
  dt$dateN[which(dt$date=="23/10/2017" | dt$date=="24/10/2017" | dt$date=="25/10/2017"
                 | dt$date=="26/10/2017" | dt$date=="27/10/2017")] <- "3 weeks"
  
  dt$dateN[which(dt$date=="04/12/2017" | dt$date=="05/12/2017")] <- "9 weeks"
  
  dt$dateN[which(dt$date=="01/03/2018" | dt$date=="05/03/2018" | dt$date=="06/03/2018"
                 | dt$date=="07/03/2018" | dt$date=="08/03/2018" | dt$date=="09/03/2018")] <- "6 months"
  
  dates <- c(na.omit(unique(dt$dateN)))
  # set as array dimension (number children, number timesteps, number max of repeats)
  children <- unique(dt$child_id)
  KK <- array(NA,dim = c(length(children),4,6))
  for (i in 1:length(children)){
    KK[i,,] <- getKKChild(children[i],dt,dates)
  }
  return(KK)
}

getKKChild <- function(ID,dt,dates){
  dts <- subset(dt,child_id==ID,select = c(child_id,Sm_A,Sm_B,dateN))
  KKChild <- matrix(NA,nrow=length(dates),ncol=6)
  for (i in 1:length(dates)){
    KKChild[i,] <- getKKChildWeek(dates[i],dts)
  }
  return(KKChild)
}

getKKChildWeek <- function(weekN,dts){
  repeatsKK <- c(dts$Sm_A[which(dts$dateN==weekN)],dts$Sm_B[which(dts$dateN==weekN)])
  length(repeatsKK)<-6 #set max number of repeats!
  return(repeatsKK)
}

getKKChildIDs <- function(nameFile){
  dt <- read.csv(nameFile)
  return(sort(unique(dt$child_id)))
}

#### Normal CCA Data ####

loadCCAdata <- function(nameFile){
  dt <- read.csv(nameFile, stringsAsFactors = F)
  dt$Poppy[which(dt$Poppy == "I*" | dt$Poppy=="-")]<-NA
  
  # restructure for categorical
  dt$Poppy[which(dt$Poppy=="T")] <- 1
  dt$Poppy[which(dt$Poppy==3)]<-4
  dt$Poppy[which(dt$Poppy==2)]<-3
  dt$Poppy[which(dt$Poppy==1)]<-2
  dt$Poppy[which(dt$Poppy==0.5)]<-1
  
  # name different weeks
  dt$dateN <- NA
  dt$dateN[which(dt$date_on_tube=="25/09/2017" | dt$date_on_tube=="26/09/2017" | dt$date_on_tube=="27/09/2017"
                 | dt$date_on_tube=="28/09/2017" | dt$date_on_tube=="29/09/2017" | dt$date_on_tube=="02/10/2017")] <- "Pre-T"
  
  dt$dateN[which(dt$date_on_tube=="23/10/2017" | dt$date_on_tube=="24/10/2017" | dt$date_on_tube=="25/10/2017"
                 | dt$date_on_tube=="26/10/2017" | dt$date_on_tube=="27/10/2017")] <- "3 weeks"
  
  dt$dateN[which(dt$date_on_tube=="04/12/2017" | dt$date_on_tube=="05/12/2017")] <- "9 weeks"
  
  dt$dateN[which(dt$date_on_tube=="01/03/2018" | dt$date_on_tube=="05/03/2018" | dt$date_on_tube=="06/03/2018"
                 | dt$date_on_tube=="07/03/2018" | dt$date_on_tube=="08/03/2018" | dt$date_on_tube=="09/03/2018")] <- "6 months"
  
  dates <- c(na.omit(unique(dt$dateN)))
  # set as array dimension (number children, number timesteps)
  children <- unique(dt$CID)
  CCA <- array(NA,dim = c(length(children),4))
  for (i in 1:length(children)){
    CCA[i,] <- getCCAChild(children[i],dt,dates)
  }
  return(CCA)
} # end of  function 

# function for getting CCA child
getCCAChild <- function(ID,dt,dates){
  dts <- subset(dt,CID==ID,select = c(CID,Poppy,dateN))
  CCAChild <- rep(NA,length(dates))
  for (i in 1:length(dates)){
    if (!is.null(which(dts$dateN==dates[i]))){
      CCAChild[i] <- ifelse(is.null(dts$Poppy[which(dts$dateN==dates[i])]),NA,as.numeric(dts$Poppy[which(dts$dateN==dates[i])]))
    }
  }
  return(CCAChild)
}

# function for getting CCA child ID
getCCAChildIDs <- function(nameFile){
  dt <- read.csv(nameFile)
  return(sort(unique(dt$CID)))
}

#### GScore data ####

loadgscoredata <- function(nameFile){
  #restructure so that it is the right scoring so that a
  # score of 1 which is negative, is given a score of 0.
  dt <- read.csv(nameFile, stringsAsFactors = F)
  
  dt$gscore <- dt$gscore-1
  
  # name different weeks
  dt$dateN <- NA
  dt$dateN[which(dt$date_on_tube=="25/09/2017" | dt$date_on_tube=="26/09/2017" | dt$date_on_tube=="27/09/2017"
                 | dt$date_on_tube=="28/09/2017" | dt$date_on_tube=="29/09/2017" | dt$date_on_tube=="02/10/2017")] <- "Pre-T"
  
  dt$dateN[which(dt$date_on_tube=="23/10/2017" | dt$date_on_tube=="24/10/2017" | dt$date_on_tube=="25/10/2017"
                 | dt$date_on_tube=="26/10/2017" | dt$date_on_tube=="27/10/2017")] <- "3 weeks"
  
  dt$dateN[which(dt$date_on_tube=="04/12/2017" | dt$date_on_tube=="05/12/2017")] <- "9 weeks"
  
  dt$dateN[which(dt$date_on_tube=="01/03/2018" | dt$date_on_tube=="05/03/2018" | dt$date_on_tube=="06/03/2018"
                 | dt$date_on_tube=="07/03/2018" | dt$date_on_tube=="08/03/2018" | dt$date_on_tube=="09/03/2018")] <- "6 months"
  
  dates <- c(na.omit(unique(dt$dateN)))
  # set as array dimension (number children, number timesteps)
  children <- unique(dt$CID)
  CCAgscore <- array(NA,dim = c(length(children),4))
  for (i in 1:length(children)){
    CCAgscore[i,] <- getCCAChildgscore(children[i],dt,dates)
  }
  return(CCAgscore)
} # end of  function 

# function for getting CCA child
getCCAChildgscore <- function(ID,dt,dates){
  dts <- subset(dt,CID==ID,select = c(CID,gscore,dateN))
  CCAChildgscore <- rep(NA,length(dates))
  for (i in 1:length(dates)){
    if (!is.null(which(dts$dateN==dates[i]))){
      CCAChildgscore[i] <- ifelse(is.null(dts$gscore[which(dts$dateN==dates[i])]),NA,as.numeric(dts$gscore[which(dts$dateN==dates[i])]))
    }
  }
  return(CCAChildgscore)
}

# function for getting CCA child ID
getCCAChildIDsgscore <- function(nameFile){
  dt <- read.csv(nameFile)
  return(sort(unique(dt$CID)))
}

#### Supplementary Material Functions ####

# this is just the function for the KK data with reinfection happeneing between BL and 3wks

do.KK.run.t2 <- function (N, Ti, R, KK) {
  ## Set seed ##
  .RNG.seed <- function(chain)
    return( switch(chain, "1"= 1, "2"= 2) )
  .RNG.name <- function(chain)
    return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )
  
  ## Initialize Status ##
  IniStatus <- rep(1,N)
  Clearance <- matrix(0,nrow=N,ncol=Ti-1)
  Reinfec <- matrix(1,nrow=N,ncol=Ti-1)
  
  m <- "model {
  
  # Prior prevalence / clearance / reinfection #
  
  prev ~ dbeta(1,1)
  
  clearance[2] ~ dbeta(1,1)
  clearance[3] ~ dbeta(1,1)
  clearance[4] ~ dbeta(1,1)
  
  reinfec[2] ~ dbeta(1,1)
  reinfec[3] ~ dbeta(1,1)
  reinfec[4] ~ dbeta(1,1)
  
  # Prior KKs #
  
  rtnb ~ dgamma(286.3235, 223.4650)
  sh ~ dgamma(83.88592,125.70331)
  rt1 ~ dbeta(47.13542,633.08366)
  rt <- rt1/(1-rt1)
  
  # accounting for the difference in length of time steps because there will be different amounts of variance
  # these are the priors for the random walk 
  sigma[4] ~ dgamma(0.001, 0.001) # sigma when the transition from timesteps up to 4
  sigma[3] <- sigma[4]*6/15 # sigma when the transition from timestep 2 up to 3 scaling down the variance from the previous step
  sigma[2] <- sigma[4]*3/15 # sigma when the transition from timestep 
  
  for(n in 1:N){	# run through pop
  
  IniStatus[n] ~ dbern(prev) # the probability of you being infected at the beginning is dependent on the prevalence of infection, doesn't this account for some of them being infected and some of them not at baseline because some zero counts might still be infected at low level infections and so might show up as a 1 or 2 on the CCA scale 
  Status[n,1] <- IniStatus[n] # individual n's status at t1 is taken from the draw above
  
  # this are the true baseline counts set outside of the time loop for the auto-regressive part 
  tKK[n,1,1] <- 0
  tKK[n,2,1] ~ dgamma(sh,rt) 
  
  for (t in 1:Ti){ # run through timepoints from timept 2 onwards
  lambda[n,t] <- tKK[n,Status[n,t]+1,t] # this needs to go here because each time the r loop reiterates it will replace this value
  for(r in 1:R){ # run through repeated measures to set the baseline KK over the repeated measures, 
  KK[n,t,r] ~ dnegbin(rtnb/(lambda[n,t]+rtnb),rtnb)
  } # end or r loop
  } # end timestep 1 t loop
  
  for (t in 2:Ti){ # HMM component
  tKK[n,1,t] <- 0
  tKK[n,2,t] ~ dnorm(tKK[n,Status[n,t-1]+1,t-1],sigma[t]) # taking the mean in the previous time step + some sort of movement around that mean, the reason it is a normal dist is because this is the movement around the previous mean
  
  Clearance[n,t-1] ~ dbern(clearance[t])
  Reinfec[n,t-1] ~ dbern(reinfec[t])
  Status[n,t] <- Status[n,t-1] * (1-Clearance[n,t-1]) + (1-Status[n,t-1]) * Reinfec[n,t-1]
  
  switch[n,t] <- abs(Status[n,t]-Status[n,t-1])      # this is a dummy variable looking at how many times an individual switches status over the run
  
  } # end hmm time
  } # end n loop
  
  #inits# .RNG.seed, .RNG.name, IniStatus, Clearance, Reinfec   
  #data# N, Ti, R, KK
  #monitor# prev, clearance, reinfec, rtnb, sh, rt, Status
}"
  
  
  # Run model #
  Results <- run.jags(m, burnin=5000, sample=10000, n.chains=2, jags.refresh = 1, method = 'parallel',
                      plots = F, silent.jags = F)
  return(Results)
}

#### do.run This second function has both CCA and KK in it ####

do.run.t2 <- function (N, Ti, R, KK, CCA) {
  ## Set seed ##
  ##this is so that they start on different values ##
  .RNG.seed <- function(chain)
    return( switch(chain, "1"= 1, "2"= 2) )
  .RNG.name <- function(chain)
    return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )
  
  ## Initialize Status ##
  IniStatus <- rep(1,N)
  Clearance <- matrix(0,nrow=N,ncol=Ti-1)
  Reinfec <- matrix(1,nrow=N,ncol=Ti-1)
  
  m <- "model {
  # Prior prevalence / clearance / reinfection #
  prev ~ dbeta(1,1)
  
  clearance[2] ~ dbeta(1,1)
  clearance[3] ~ dbeta(1,1)
  clearance[4] ~ dbeta(1,1)
  
  
  reinfec[2] ~ dbeta(1,1)
  reinfec[3] ~ dbeta(1,1)
  reinfec[4] ~ dbeta(1,1)
  
  # Prior KKs 
  rtnb ~ dgamma(286.3235, 223.4650)
  sh ~ dgamma(0.001,0.001)
  rt1 ~ dbeta(1,1)
  rt <- rt1/(1-rt1)
  
  # Priors for the random walk #
  # accounting for the difference in length of time steps because there will be different amounts of variance
  sigma[4] ~ dgamma(0.001, 0.001) # sigma when the transition from timesteps up to 4
  sigma[3] <- sigma[4]*6/15 # sigma when the transition from timestep 2 up to 3 scaling down the variance from the previous step
  sigma[2] <- sigma[4]*3/15 # sigma when the transition from timestep 
  
  # Priors CCA #
  k ~ dnorm(0.063291330, 1/0.009817479^2)
  intercept ~ dunif(0.0139464, 8.5045100)
  TracPosProb ~ dbeta(4.379053,25.081881)
  
  # MODEL 
  for(n in 1:N){	# run through pop
  
  IniStatus[n] ~ dbern(prev)
  Status[n,1] <- IniStatus[n]
  
  # these are the true baseline counts set outside of the time loop for the auto-regressive part 
  tKK[n,1,1] <- 0
  tKK[n,2,1] ~ dgamma(sh,rt) 
  
  prob[n,1,1] ~ dbern(TracPosProb)
  prob[n,2,1] <- 1 / (1 + exp(-(tKK[n,2,1])))
  
  CCA[n,1] ~ dbinom(prob[n,Status[n,1]+1,1],(3*Status[n,1]+1))
  
  for (t in 1:Ti){ # run through timepoints
  lambda[n,t] <- tKK[n,Status[n,t]+1,t] # this needs to go here because each time the r loop reiterates it will replace this value
  for( r in 1:R){  # run through repeat measures
  KK[n,t,r] ~ dnegbin(rtnb/(lambda[n,t]+rtnb),rtnb) # generating the data with noise and then sampling from the dataset with the gamma sh/rt parameters?
  } # end or r loop
  } # end timestep 1 t loop
  
  for (t in 2:Ti){ # HMM component
  tKK[n,1,t] <- 0                                                                                                                                                                                                                                                                                               
  tKK[n,2,t] ~ dnorm(tKK[n,Status[n,t-1]+1,t-1],sigma[t]) 
  
  prob[n,1,t] ~ dbern(TracPosProb)
  prob[n,2,t] <- 1 / (1 + exp(-k*(tKK[n,2,t]-intercept))) # logistic function, for a given KK
  # probabliliy that feeds into the binomial for the CCA below, so that you have the 4 draws of whether you are positive or not. 
  # it includes the trueKK because the logistic function describes the point at which the KK test can pick up infection dependent the worm burden as determined from the KK count scaled by k) 
  
  CCA[n,t] ~ dbinom(prob[n,Status[n,t]+1,t],(3*Status[n,t]+1))
  
  Clearance[n,t-1] ~ dbern(clearance[t])
  Reinfec[n,t-1] ~ dbern(reinfec[t])
  
  Status[n,t] <- Status[n,t-1] * (1-Clearance[n,t-1]) + (1-Status[n,t-1]) * Reinfec[n,t-1] 
  
  switch[n,t] <- abs(Status[n,t]-Status[n,t-1])      # this is a dummy variable looking at how many times an individual switches status over the run
  }
  
  }  
  
  #inits# .RNG.seed, .RNG.name, IniStatus, Clearance, Reinfec         
  #data# N, Ti, R, KK, CCA
  #monitor#  prev, clearance, reinfec, rtnb, TracPosProb, k, intercept, sh, rt, Status
}"
  
  # Run model #
  Results <- run.jags(m, burnin=5000, sample=10000, n.chains=2, jags.refresh = 1, method = 'parallel',
                      plots = F, silent.jags = F)
  return(Results)
}

#### Model Output Functions ####

#### make whole dataset into factors ####

df.factor <- function(x){
  x[sapply(x, is.numeric)] <- lapply(x[sapply(x, is.numeric)], as.factor)  
}

#### schisto status manipulation ####
# this is to look at the proportion of time each individual spends infected in each time step. 

time.steps <- function(model.output){

  t1 <- as.data.frame(model.output[,1:210])
  t2 <- as.data.frame(model.output[,211:420])
  t3 <- as.data.frame(model.output[,421:630])
  t4 <- as.data.frame(model.output[,631:840])
  
  t1[nrow(t1)+1,] <- colSums(t1)
  t1[nrow(t1),] <- t1[nrow(t1),]/40000
  
  t1.means <- as.data.frame(t1[nrow(t1),1:210])
  
  t1.means$time <-as.factor("Baseline")
  
  T1 <- melt(t1.means)
  
  t2[nrow(t2)+1,] <- colSums(t2)
  t2[nrow(t2),] <- t2[nrow(t2),]/40000
  
  t2.means <- as.data.frame(t2[nrow(t2),1:210])
  
  t2.means$time <-as.factor("ThreeWeeks")
  
  T2 <- melt(t2.means)
  
  t3[nrow(t3)+1,] <- colSums(t3)
  t3[nrow(t3),] <- t3[nrow(t3),]/40000
  
  t3.means <- as.data.frame(t3[nrow(t3),1:210])
  
  t3.means$time <- as.factor("NineWeeks")
  
  T3 <- melt(t3.means)
  
  t4[nrow(t4)+1,] <- colSums(t4)
  t4[nrow(t4),] <- t4[nrow(t4),]/40000
  
  t4.means <- as.data.frame(t4[nrow(t4),1:210])
  
  t4.means$time <- as.factor("SixMonths")
  
  T4 <- melt(t4.means)
  
  props <- rbind(T1, T2, T3, T4)
  
  props$time <- as.factor(props$time)
  props$variable <- as.factor(props$variable)
  
  return(props)
}

#### standard error and 95% confidence interval function ####

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  #datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#### add legend fuction ####
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             # Original function by Jan Gleixner (@jan-glx)
                             # Adjustments by Wouter van der Bijl (@Axeman)
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                               quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           }
)

create_quantile_segment_frame <- function(data, draw_quantiles, split = FALSE, grp = NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(
      x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  } else {
    data.frame(
      x = ggplot2:::interleave(violin.xminvs, violin.xs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  }
}

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

# putting letter labels on the facet figures 

tag_facet2 <- function(p, open = "", close = "", tag_pool = LETTERS, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}
