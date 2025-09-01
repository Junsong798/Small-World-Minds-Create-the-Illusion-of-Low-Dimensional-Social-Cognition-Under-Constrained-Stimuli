# Trait Covariance Evolution Simulation
# BIAS 3 model


######Packages######
library(psych)
library(plyr)






######Parameters######

#Model Loops#
#Set the total number of times the model will run
modelloops<-100



#Population Size#
#Set the population size
popsize<-200






######Functions######


#Agent generation#

agentgenerate<-function(popsize,selfphysatt,matephysatt,sex){
  if(popsize>0){
    
    #Generate "perceptions" of self and mate physical attractiveness
    #These perceptions are reasonably accurate--actual trait value explains 60% of the variance in trait perception
    selfphysatt<-selfphysatt+rnorm(length(selfphysatt),0,sqrt((sd(selfphysatt)^2-.60*sd(selfphysatt)^2)/.60))
    matephysatt<-matephysatt+rnorm(length(matephysatt),0,sqrt((sd(matephysatt)^2-.60*sd(matephysatt)^2)/.60))
    
    #Physical attractiveness perceptions become the source of bias for all other traits
    #Determine the amount of noise to add to physical attractiveness perception to create the biased perception of traits
    #This bias will account for approximately 45% of the variance in traits overall by default (close to what is observed in the primary model)
    selfnoisesd<-sqrt((sd(selfphysatt)^2-.45*sd(selfphysatt)^2)/.45)
    matenoisesd<-sqrt((sd(matephysatt)^2-.45*sd(matephysatt)^2)/.45)
    
    selfnoisesd<-seq(.05,selfnoisesd*2,length.out=9)
    matenoisesd<-seq(.05,matenoisesd*2,length.out=9)
    
    
    #Generate self-trait perceptions by adding the appropriate noise to physical attractiveness perceptions
    traits<-cbind(selfphysatt,sapply(1:9,function(x) selfphysatt+rnorm(popsize,0,selfnoisesd[x])))
    colnames(traits)<-paste(paste('trait',1:ncol(traits),sep=""))
    
    #Agent preferences are set to 7 for all traits
    #This simplification ensures that bias manipulates traits directly with respect to desirability
    preferences<-matrix(7,popsize,10)
    colnames(preferences)<-paste(paste('preference',1:ncol(preferences),sep=""))
    
    #Generate mate trait perceptions for each agent by adding noise to the physical attractiveness perception
    mate<-cbind(matephysatt,sapply(1:9,function(x) matephysatt+rnorm(popsize,0,matenoisesd[x])))
    colnames(mate)<-paste(paste('mtrait',1:ncol(mate),sep=""))
    
    sex<-sex
    
    return(as.data.frame(cbind(traits,preferences,mate,sex)))
  }
}







######Model Start######

#Loop through the preset number of model loops
for(m in 1:modelloops){
  
  #Generate male physical attractiveness
  mphysatt<-rnorm(popsize/2,4,2)
  
  #Generate female physical attractiveness so as to be correlated with male attractiveness
  fphysatt<-mphysatt+rnorm(popsize/2,0,sqrt((sd(mphysatt)^2-.09*sd(mphysatt)^2)/.09))
  
  #Rescale female physical attractiveness to the same common scale
  fphysatt<-4+scale(fphysatt)*2
  
  #Generate agents
  males<-agentgenerate(popsize/2,mphysatt,fphysatt,1)
  females<-agentgenerate(popsize/2,fphysatt,mphysatt,0)
  
  #Give each agent a unique PIN
  males$PIN<-sample(1:nrow(males),nrow(males),replace=F)
  females$PIN<-sample(nrow(males):(nrow(males)+nrow(females)),nrow(females),replace=F)
  
  datam<-males
  dataf<-females
  
  data<-rbind(dataf,datam)
  
  
  #Create a unique filename so I don't overwrite old results.
  #Name format is "Algorithm Frequencies MonthDayYear HourMinute"
  path<-"C:/Users/conroy-beam/Google Drive/Research/Trait Covariance Evolution/Simulation Outputs/BIAS3/Couple Population/TCE Couple Population M"
  
  format<-".csv"
  date<-format(Sys.time(),format="%m%d%Y %H%M%S")
  file<-file.path(paste0(path,m," ",date,format))
  
  write.csv(data,file=file)
  
}
