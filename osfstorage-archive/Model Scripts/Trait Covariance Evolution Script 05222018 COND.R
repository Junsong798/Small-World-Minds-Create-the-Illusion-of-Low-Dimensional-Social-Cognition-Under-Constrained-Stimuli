# Trait Covariance Evolution Simulation
# 
# ###Purpose###
# An agent-based model exploring the evolution of trait covariation in response to assortative mating. In this model, traits develop based on an inherited condition value and mate choice is random.
# 
# ###Agent description###
# 200 agents will be generated. Each agent will have 10 mate preferences and corresponding traits, a sex, and a unique identifying number (PIN).
# 
# Traits and preferences
# Each agent will have 10 traits. Traits will be generated from random normal distributions.
# Agents will also have preferences. Preferences will be generated randomly from distributions identical to those used to generate traits.
# 
# Sex
# All agents will be randomly assigned to be either male or female.
#
# ###Life cycle###
# 1. Compute attraction
# 2. Select mates
# 3. Reproduce
# 4. Die
# 
# Computing attraction
# Each agent will compute their attraction to all opposite-sex agents as the Euclidean distance between the agent's own preferences and the traits of each potential mate.
# These distances will be scaled and transformed such that a value of 10 means a perfect match between the agent's preferences and the potential mate's traits and a value of 0 means worst possible match.
# 
# Mate selection
# The simulation will next compute the mutual attraction matrix for all agents. This is the product of attraction values for all possible couples. The selection procedure will next follow several steps. 
# (1) The model will select a random couple from the mutual attraction matrix.
# (2) These two agents will be paired and placed into a separate matrix of paired agents. 
# (3) The selection procedure will remove the paired male and female from the mutual attraction matrix.
# This process will iterate until all possible pairs are made. If there are agents who cannot pair with a mate (e.g. unequal sex ratios), these agents will not be placed in the paired agent matrix.
# 
# Reproduction
# Paired couples will reproduce. Couples will produce offspring in proportion to the sum of their "energy" values, where each agent's energy is proportional to the distance between each of their trait values and the optimum trait value.
# Offspring will inherit traits and preferences from their parents. For each preference, offspring will randomly inherit the value of either their mother or father. Offspring will inherit their condition value randomly from either of their parents. This condition value will in turn determine each of their traits.
# A small amount of random noise will next be added to each inherited condition and preference value to simulate mutation.
# 
# Death
# Paired couples will reproduce. Couples will produce offspring in proportion to the sum of their "energy" values, where each agent's energy is proportional to the distance between each of their trait values and pre-set optimum trait values.
# At this stage, the full life cycle is complete. The next generation will start at the beginning of the life cycle.
# This will repeat for 1000 generations of evolution


######Packages######
library(psych)
library(plyr)
library(foreach)
library(doParallel)




######Parameters######

#Model Loops#
#Set the total number of times the model will run
modelloops<-100



#Population Size#
#Set the population size
popsize<-200



#Mutation Rate#
#Set the variability added to traits and preferences each generation
mr<-.006



#Generations#
#Number of generations of evolution
generations<-1000



#Strength of Selection#
#Set the strength of selection
#This value will be proportion of variance in probability of reproducing that is explained by couples' pooled energy values
selstrength<-.10


#Number of Cores#
#Number of cores to use for the parallel loop
clust<-makeCluster(detectCores()-1)
registerDoParallel(clust)






######Preset File Locations######

#Factor Variance Data#
#Create a filepath for saving the factor variance data
#Create a unique filename so I don't overwrite old results.
#Name format is "TCE Factor Variance MonthDayYear HourMinute"
fvpath<-choose.dir(caption="Select Folder for Factor Variance Data")
fvpath<-paste0(fvpath,"\\TCE COND Factor Variance ")

#Couple Population Data#
#Create a filepath for saving the couple population data

#Create a unique filename so I don't overwrite old results.
#Name format is "TCE Couple Population MonthDayYear HourMinute"
cppath<-choose.dir(caption="Select Folder for Couple Population Data")
cppath<-paste0(cppath,"\\TCE COND Couple Population ")








######Functions######


#Agent generation#

agentgenerate<-function(popsize,sex,fitvec){
  if(popsize>0){
    
    #Each agent is assigned a "condition" variable
    condition<-rnorm(popsize,4,2)
    
    #Determine the amount of noise to add to the condition variable to create traits
    #Condition will explain roughly 45% of the variance in traits overall by default (close to what is observed in the primary model)
    noisesd<-sqrt((sd(condition)^2-.45*sd(condition)^2)/.45)
    
    #Generate traits by adding the appropriate noise to condition
    traits<-apply(matrix(,popsize,10),2,function(x) condition+rnorm(popsize,0,noisesd))
    colnames(traits)<-paste(paste('trait',1:ncol(traits),sep=""))
    
    #Agent preferences are drawn from random normal distributions. Preferences are named appropriately and contrained to values between 1 and 7
    preferences<-apply(matrix(,popsize,10),2,function(x) rnorm(popsize,4,2))
    colnames(preferences)<-paste(paste('preference',1:ncol(preferences),sep=""))
    
    #Calculate energy for all agents
    #Computes a weight for each optimum value in the fitness vector
    #This weight adjusts for variability in the range of possible deviations across optimum values
    #Without this weight, selection doesn't act as strongly on traits where the optimum value falls toward the center of the trait range
    fitvecweights<-sapply(fitvec,simplify=T,function(x) 10/max(abs(x-1:7)))
    
    #Calculate agent energy based on deviation from optimum using the above weights
    energy<-apply(traits,1,function(x) sum(10-(abs(fitvec-x)*fitvecweights)))
    
    sex<-sex
    
    return(as.data.frame(cbind(traits,preferences,condition,energy,sex)))
  }
}



#Attraction calculation#
attraction<-function(agentprefs,targettraits){
  #agentprefs is the ideal preferences of an individual agent
  #targettraits is the traits of all opposite-sex agents
  
  #Renames the agent preferences and traits so they can be rbound
  agentprefs<-agentprefs[1:10]
  names(agentprefs)<-1:10
  colnames(targettraits)<-1:10
  
  #Calculates the Euclidean distance between the agent's preferences and each mate's traits
  #This is caled and transformed so that a value of 10 means shortest possible distance and a value of 0 means longest possible distance
  att<-apply(targettraits,1,function(x) 10*((-1*dist(rbind(agentprefs,x)))+sqrt(6^2*10))/(sqrt(6^2*10)))
  
  return(att)
  
}


#Reproduction#
reproduce<-function(repro,mr,fitvec){
  
  #Empty matrix for storing offspring
  offspring<-matrix(NA,nrow(repro),23)
  
  #Iterate through each couple, producing an offspring for each iteration
  for(o in 1:nrow(repro)){
    #Separate out the current mother and father
    mother<-repro[o,1:24]
    father<-repro[o,25:48]
    
    #Offspring inherit their condition value randomly from either parent
    condition<-rbinom(1,1,.5)
    condition<-ifelse(condition==1,as.numeric(father)[21],as.numeric(mother)[21])
    
    #Add a small amount of variability to condition
    condition<-condition+rnorm(1,0,.06)
  
    #Temporarily create trait values as equivalent to condition values
    traits<-rep(condition,10)
    
    #Have agents inherit preferences
    preferences<-rbinom(10,1,.5)
    preferences<-ifelse(preferences==1,as.numeric(father[11:20]),as.numeric(mother[11:20]))
    
    #Add a set amount of variability to preferences
    preferences<-preferences+rnorm(length(preferences),0,mr)
    
    #Create the offspring, with temporary NA values for sex and energy
    offspring[o,]<-c(traits,preferences,condition,NA,NA)
  }
  
  #Determine the amount of noise to add to the condition variable to create traits
  #Condition will explain roughly 80% of the variance in traits by default (close to what is observed in the primary model)
  noisesd<-sqrt((sd(offspring[,21])^2-.45*sd(offspring[,21])^2)/.45)
  
  #Add noise to traits
  traits<-apply(offspring[,1:10],2,function(x) x+rnorm(nrow(offspring),0,noisesd))
  
  #Save the noisy traits over the offspring's condition values
  offspring[,1:10]<-traits
  
  #Calculate the offspring's actual energy
  #Computes a weight for each optimum value in the fitness vector
  #This weight adjusts for variability in the range of possible deviations across optimum values
  #Without this weight, selection doesn't act as strongly on traits where the optimum value falls toward the center of the trait range
  fitvecweights<-sapply(fitvec,simplify=T,function(x) 10/max(abs(x-1:7)))
  
  for(o in 1:nrow(offspring)){
    
    #Calculate offspring energy based on deviation from optimum using the above weights
    offspring[o,22]<-sum(10-(abs(fitvec-offspring[o,1:10])*fitvecweights))
  }
  
  #Convert the offspring matrix to a dataframe
  offspring<-as.data.frame(offspring)
  
  #Name everything for ease of use
  colnames(offspring)<-colnames(repro)[1:23]
  
  #Assign each offspring a random sex
  offspring$sex<-sample(rep(c(0,1),each=(nrow(offspring)/2)))
  
  return(offspring)
}












######Model Start######

#Loop through the preset number of model loops
foreach(m=1:modelloops,.packages=c("plyr","psych")) %dopar% {
  
  ######Agent Generation######
  
  #Generate the fitness vector
  #This vector determines what trait value will return the agent the maximum amount of energy
  #Initially sets the optimum trait value to 7 for all traits, but this can be set to whatever.
  fitvec<-rep(7,10)
  
  #Generate agents
  males<-agentgenerate(popsize/2,1,fitvec)
  females<-agentgenerate(popsize/2,0,fitvec)
  
  #Give each agent a unique PIN
  males$PIN<-sample(1:nrow(males),nrow(males),replace=F)
  females$PIN<-sample(nrow(males):(nrow(males)+nrow(females)),nrow(females),replace=F)
  
  
  ########## Life Cycle ##########
  #Life cycle begins here
  #Iterations at this level represent generations, set above.
  
  #Create a dataframe for storing variance explained by the d-factor
  facvar<-matrix(,1,(generations+1))
  
  datam<-males
  dataf<-females
  
  #calculates the average preferences of each sex
  femaleprefs<-colMeans(dataf[,11:20])
  maleprefs<-colMeans(datam[,11:20])
  
  #Converts trait values on each dimension to trait desirabilities on each dimension
  datam2<-data.frame(t(apply(datam[,1:10],1,function(x) 6-abs(x-femaleprefs))))
  dataf2<-data.frame(t(apply(dataf[,1:10],1,function(x) 6-abs(x-maleprefs))))
  
  #Extracts a general desirability factor from desirabilities along each dimension
  mfa<-principal(datam2,nfactors=1)
  ffa<-principal(dataf2,nfactors=1)
  
  #Calculate the proportion of variance accounted for by the general desirability factor
  mvar<-colSums(mfa$loading*mfa$loading)/dim(mfa$loading)[1]
  fvar<-colSums(ffa$loading*ffa$loading)/dim(ffa$loading)[1]
  
  facvar[1,1]<-mean(c(mvar,fvar))
  
  
  
  for(g in 1:generations){
    
    
    
    ######Computing Attraction######
    
    #Males
    #Calculate how attracted each male is to each female
    maleattmatrix<-t(apply(males,1,function(q) attraction(q[11:20],females[,1:10])))
    
    
    #Females
    #Calculate how attracted each female is to each male
    femaleattmatrix<-t(apply(females,1,function(q) attraction(q[11:20],males[,1:10])))
    
    
    
    ######Mate Selection######
    
    #Mutual attraction matrix
    #Computes the mutual attraction matrix--how mutually attracted each possible couple would be.
    mutattmatrix<-as.data.frame(maleattmatrix*t(femaleattmatrix))
    
    #Renames the columns to be the female PIN values
    colnames(mutattmatrix)<-females$PIN
    
    #Renames the rows to be male PIN values
    rownames(mutattmatrix)<-males$PIN
    
    #A blank dataframe for storing the paired couples
    pairs<-matrix(,min(nrow(males),nrow(females)),(ncol(males)+ncol(females)))
    pairs<-as.data.frame(pairs)
    
    #Iterates through the mutual attraction matrix
    for(p in 1:(popsize/2)){
      
      #Determines all locations in the mutual attraction matrix in which neither agent has been paired
      maxmut<-which(is.na(mutattmatrix)==F,arr.ind=T)
      
      #If there is more than one valid location, picks a random location
      if(length(maxmut)/2 > 1){maxmut<-maxmut[sample(nrow(maxmut),1,replace=F),]}
      
      #Determines which male and which female were actually paired
      pairedmale<-males[males$PIN==rownames(mutattmatrix)[maxmut[1]],]
      pairedfemale<-females[females$PIN==colnames(mutattmatrix)[maxmut[2]],]
      
      #Places the paired male and female in the pairs dataframe
      pairs[p,]<-cbind(pairedfemale,pairedmale)
      
      #Removes the paired male and female from the mutual attraction matrix
      mutattmatrix[maxmut[1],]<-NA
      mutattmatrix[,maxmut[2]]<-NA
      
    }
    
    #Just names the pairs matrix appropriately so it's easier to read. Probably unnecessary.
    colnames(pairs)<-c(colnames(females),colnames(males))
    
    ######Reproduce######

    #Determine each couple's summed energy values
    pairse<-(pairs[,22]+pairs[,46])
    
    #Rescale these values just to make imposing selection strength easier
    pairse<-pairse-min(pairse)
    pairse<-100*pairse/max(pairse)
    
    #Add the selstrength constant, which controls the strength of selection on energy values
    pairse<-pairse+(100/selstrength)
    
    #Determines how often each couple gets to reproduce by sampling randomly from couples in proportion to their energy values
    repro<-pairs[sample(nrow(pairs),popsize,replace=T,prob=pairse),]
    
    #Have all agent couples produce offspring
    offspring<-reproduce(repro,mr,fitvec)
    
    ######Death######
    #Give each offspring a PIN
    offspring$PIN<-sample(1:nrow(offspring),nrow(offspring),replace=F)
    
    #Save the offspring over the parent dataframes. This kills the parents
    males<-subset(offspring,offspring$sex==1)
    females<-subset(offspring,offspring$sex==0)
    
    ######Data Saving######
    #Here is where I save whatever data is of interest
    
    #Split the last pairs dataframe into males and females
    datam<-pairs[,25:48]
    dataf<-pairs[,1:24]
    
    #First pulls out males' preferences.
    datampref<-datam[,11:20]
    
    #Then pulls males' traits
    datamself<-datam[,1:10]
    
    #Pulls out females' preferences
    datafpref<-dataf[,11:20]
    
    #Then pulls out females' traits 
    datafself<-dataf[,1:10]
    
    #Renames the matrices so they can be rbound
    names(datampref)<-1:10
    names(datamself)<-1:10
    
    names(datafpref)<-1:10
    names(datafself)<-1:10
    
    #calculates the average preferences of each sex
    femaleprefs<-colMeans(datafpref)
    maleprefs<-colMeans(datampref)
    
    #Calculates the mate value of each male. The Euclidean distance between their traits and the preferences of the opposite sex
    datam$selfmv<-apply(datamself,1,function(x) dist(rbind(x,femaleprefs)))
    datam$selfmv<-10*((-1*datam$selfmv)+(sqrt(6^2*10)))/(sqrt(6^2*10))
    
    #Calculates the mate value of each female. The Euclidean distance between their traits and the preferences of the opposite sex
    dataf$selfmv<-apply(datafself,1,function(x) dist(rbind(x,maleprefs)))
    dataf$selfmv<-10*((-1*dataf$selfmv)+(sqrt(6^2*10)))/(sqrt(6^2*10))
    
    #Converts trait values on each dimension to trait desirabilities on each dimension
    datam2<-data.frame(t(apply(datam[,1:10],1,function(x) 6-abs(x-femaleprefs))))
    dataf2<-data.frame(t(apply(dataf[,1:10],1,function(x) 6-abs(x-maleprefs))))
    
    #Extracts a general desirability factor from desirabilities along each dimension
    mfa<-principal(datam2,nfactors=1)
    ffa<-principal(dataf2,nfactors=1)
    
    #Calculate the proportion of variance accounted for the general desirability factor
    mvar<-colSums(mfa$loading*mfa$loading)/dim(mfa$loading)[1]
    fvar<-colSums(ffa$loading*ffa$loading)/dim(ffa$loading)[1]
    
    facvar[1,(g+1)]<-mean(c(mvar,fvar))
    
    datam$fa<-as.numeric(mfa$scores)
    dataf$fa<-as.numeric(ffa$scores)
    
    data<-rbind(dataf,datam)
    
  }
  

  #Save the factor variance data
  format<-".csv"
  date<-format(Sys.time(),format="%m%d%Y %H%M%S")
  fvfile<-file.path(paste0(fvpath,"M",m," ",date,format))
  
  write.csv(facvar,file=fvfile)
  
  
  #Save the couple population data
  format<-".csv"
  date<-format(Sys.time(),format="%m%d%Y %H%M%S")
  cpfile<-file.path(paste0(cppath,"M",m," ",date,format))
  
  write.csv(data,file=cpfile)
}



stopCluster(clust)


