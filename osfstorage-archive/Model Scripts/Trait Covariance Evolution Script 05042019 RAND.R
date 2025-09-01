# Trait Covariance Evolution Simulation
# 
# ###Purpose###
# An agent-based model exploring the evolution of trait covariation in response to assortative mating. In this model, traits are inherited directly and mate choice is preference-driven and assortative.
# 
# ###Agent description###
# 200 agents will be generated. Each agent will have 10 mate preferences and corresponding traits, a sex, and a unique identifying number (PIN).
# 
# Traits and preferences
# Each agent will have 10 traits. Traits will be the sum of inherited "gene" values, themselves drawn from random uniform distributions.
# Agents will also have preferences. Preferences, like traits, will be the sum of gene values drawn from random uniform distributions.
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
# (1) The most mutually attracted couple will first be identified. If more than one couple is equally mutually attracted, a random couple will be chosen.
# (2) These two agents will be paired and placed into a separate matrix of paired agents. 
# (3) The selection procedure will remove the paired male and female from the mutual attraction matrix.
# This process will iterate until all possible pairs are made. If there are agents who cannot pair with a mate (e.g. unequal sex ratios), these agents will not be placed in the paired agent matrix.
# 
# Reproduction
# Paired couples will reproduce. Couples will produce offspring in proportion to the sum of their "energy" values, where each agent's energy is proportional to the distance between each of their trait values and pre-set optimum trait values.
# Offspring will inherit trait and preference genes from their parents. For each trait and preference, offspring will randomly inherit the value of either their mother or father.
# A small amount of random noise will next be added to each inherited trait and preference value to simulate mutation.
# 
# Death
# After all agent couples reproduce, the offspring matrices will be saved over the parent matrices. The couple matrix will also be emptied. This effectively kills all parents. 
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
mr<-(.006/100)



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
fvpath<-paste0(fvpath,"\\TCE RAND Factor Variance ")

#Couple Population Data#
#Create a filepath for saving the couple population data

#Create a unique filename so I don't overwrite old results.
#Name format is "TCE Couple Population MonthDayYear HourMinute"
cppath<-choose.dir(caption="Select Folder for Couple Population Data")
cppath<-paste0(cppath,"\\TCE RAND Couple Population ")





######Functions######

#Gene Generation#

genegenerate<-function(popsize){
  if(popsize>0){
    #Trait genes are drawn from random uniform distributions and named appropriately.
    traits<-apply(matrix(,popsize,1000),2,function(x) runif(popsize,0,0.08))
    colnames(traits)<-paste(paste('trait_gene',1:ncol(traits),sep=""))
    
    #Preference genes are drawn from random uniform distributions and named appropriately.
    preferences<-apply(matrix(,popsize,1000),2,function(x) runif(popsize,0,0.08))
    colnames(preferences)<-paste(paste('preference_gene',1:ncol(preferences),sep=""))
    
    return(as.data.frame(cbind(traits,preferences)))
  }
}



#Agent generation#

agentgenerate<-function(agentgenes,sex,fitvec){
  
  #Split genes into trait and preference genes for ease of use
  traitgenes<-agentgenes[,1:1000]
  prefgenes<-agentgenes[,1001:2000]
  
  #Determine each agent's traits by summing their relevant genes and then name them appropriately.
  traits<-sapply(1:10,function(x) rowSums(traitgenes[,(1+(x-1)*100):(x*100)]))
  colnames(traits)<-paste(paste('trait',1:ncol(traits),sep=""))
  
  #Do the same for preferences
  preferences<-sapply(1:10,function(x) rowSums(prefgenes[,(1+(x-1)*100):(x*100)]))
  colnames(preferences)<-paste(paste('preference',1:ncol(preferences),sep=""))
  
  #Calculate energy for all agents
  #Computes a weight for each optimum value in the fitness vector
  #This weight adjusts for variability in the range of possible deviations across optimum values
  #Without this weight, selection doesn't act as strongly on traits where the optimum value falls toward the center of the trait range
  fitvecweights<-sapply(fitvec,simplify=T,function(x) 10/max(abs(x-1:7)))
  
  #Calculate agent energy based on deviation from optimum using the above weights
  energy<-apply(traits,1,function(x) sum(10-(abs(fitvec-x)*fitvecweights)))
  
  sex<-sex
  
  return(as.data.frame(cbind(traits,preferences,energy,sex)))
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


#Gene Reproduction#
reproducegenes<-function(mother,father,mr){
  
  #Offspring inherit trait and preference genes randomly from either of their parents
  traits<-rbinom(1000,1,.5)
  traits<-ifelse(traits==1,as.numeric(father),as.numeric(mother))
  
  preferences<-rbinom(1000,1,.5)
  preferences<-ifelse(preferences==1,as.numeric(father[1001:2000]),as.numeric(mother[1001:2000]))
  
  #Add a set amount of variability to genes
  traits<-traits+rnorm(length(traits),0,mr)
  preferences<-preferences+rnorm(length(preferences),0,mr)
  
  return(c(traits,preferences))
}


#Offspring Reproduction#
reproduceoffspring<-function(offgenes,fitvec){
  
  #Empty matrix for storing offspring
  offspring<-matrix(NA,nrow(repro),22)
  
  #Split genes in to trait and preference genes for ease of use
  traitgenes<-offgenes[,1:1000]
  prefgenes<-offgenes[,1001:2000]
  
  #Determine each offspring's traits by summing the relevant genes
  traits<-sapply(1:10,function(x) rowSums(traitgenes[,(1+(x-1)*100):(x*100)]))
  
  #Do the same for preferences
  preferences<-sapply(1:10,function(x) rowSums(prefgenes[,(1+(x-1)*100):(x*100)]))
  
  #Create the offspring, with temporary NA values for sex and energy
  offspring[,1:20]<-cbind(traits,preferences)
  
  #Calculate the offspring's actual energy
  #Computes a weight for each optimum value in the fitness vector
  #This weight adjusts for variability in the range of possible deviations across optimum values
  #Without this weight, selection doesn't act as strongly on traits where the optimum value falls toward the center of the trait range
  fitvecweights<-sapply(fitvec,simplify=T,function(x) 10/max(abs(x-1:7)))
  
  #Calculate offspring energy based on deviation from optimum using the above weights
  offspring[,21]<-sapply(1:nrow(offspring),function(x) sum(10-(abs(fitvec-offspring[x,1:10])*fitvecweights)))
  
  #Convert the offspring matrix to a dataframe
  offspring<-as.data.frame(offspring)
  
  #Assign each offspring a random sex
  offspring[,22]<-sample(rep(c(0,1),each=(nrow(offspring)/2)))
  
  return(offspring)
}











######Model Start######

#Loop through the preset number of model loops
foreach(m=1:modelloops,.packages=c("plyr","psych")) %dopar% {
  
  ######Agent Generation######
  
  #Generate the fitness vector
  #This vector determines what trait value will return the agent the maximum amount of energy
  fitvec<-runif(10,1,7)
  
  #Generate agent genes
  malegenes<-genegenerate(popsize/2)
  femalegenes<-genegenerate(popsize/2)
  
  #Generate agents
  males<-agentgenerate(malegenes,1,fitvec)
  females<-agentgenerate(femalegenes,0,fitvec)
  
  #Give each agent a unique PIN
  males$PIN<-sample(1:nrow(males),nrow(males),replace=F)
  females$PIN<-sample(nrow(males):(nrow(males)+nrow(females)),nrow(females),replace=F)
  
  #Store that PIN in the gene dataframes as well
  malegenes$PIN<-males$PIN
  femalegenes$PIN<-females$PIN
  
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
      
      #Determines all possible couples in the matrix
      maxmut<-which(is.na(mutattmatrix)==F,arr.ind=T)
      
      #If there is more than one possible couple, picks a random location
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
    pairse<-(pairs[,21]+pairs[,44])
    
    #Rescale these values just to make imposing selection strength easier
    pairse<-pairse-min(pairse)
    pairse<-100*pairse/max(pairse)
    
    #Add the selstrength constant, which controls the strength of selection on energy values
    pairse<-pairse+(100/selstrength)
    
    #Determines how often each couple gets to reproduce by sampling randomly from couples in proportion to their energy values
    repro<-pairs[sample(nrow(pairs),popsize,replace=T,prob=pairse),]
    
    #Generate offspring genes based on the repro dataframe
    offgenes<-t(apply(repro,1,function(x) reproducegenes(femalegenes[femalegenes$PIN==x[23],],malegenes[malegenes$PIN==x[46],],mr)))
    
    #Convert genes to a dataframe and label for ease of use later
    offgenes<-as.data.frame(offgenes)
    colnames(offgenes)<-colnames(femalegenes)[1:2000]
    
    #Produce offspring agents based on these genes
    offspring<-reproduceoffspring(offgenes,fitvec)
    
    #Rename everything for ease of use
    colnames(offspring)<-colnames(repro)[1:22]
    
    ######Death######
    #Give each offspring a PIN
    offspring$PIN<-sample(1:nrow(offspring),nrow(offspring),replace=F)
    
    #Save the offspring over the parent dataframes. This kills the parents
    males<-subset(offspring,offspring$sex==1)
    females<-subset(offspring,offspring$sex==0)
    
    #Save the genes of the offspring
    malegenes<-offgenes[offspring$sex==1,]
    femalegenes<-offgenes[offspring$sex==0,]
    
    #Store the agent PINs in the gene dataframes
    malegenes$PIN<-males$PIN
    femalegenes$PIN<-females$PIN
    
    ######Data Saving######
    #Here is where I save whatever data is of interest
    
    #Split the last pairs dataframe into males and females
    datam<-pairs[,24:46]
    dataf<-pairs[,1:23]
    
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


