######Trait Covariance Evolution Model Analysis Script######



######Packages######
library(ggplot2)
library(plyr)
library(gridExtra)
library(psych)
library(reshape2)
library(optimx)
library(lmerTest)



######Functions######
stdcoef <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}



######Parameters######

#Seed#
#Set a seed to make the stochastic analyses reproducible
set.seed(06082018)



#Maximum Euclidean Distance#
#A constant used throughout to scale Euclidean distances
maxeuc<-sqrt(6^2*10)



#Save Figures#
sf<-1

if(sf==1){
  figdir<-choose.dir()
}





###HERIT MODELS###

##HERIT AS (h=.5)##

#Factor Evolution#

#Load factor size data
setwd(choose.dir(caption="Navigate to folder containing HERIT AS (h=.5) factor variance data"))

#Collapse facdata into a single list
data<-lapply(list.files(),read.csv)

#Tag each model run with a unique modelrun number
for(i in 1:length(data)){
  data[[i]]$modelrun<-i
}

#Collapse all data into a single dataframe
sim<-NULL
sim<-do.call(rbind.fill,data)

heritas50_facevo<-data.frame("model"=rep("Assortative",1001),"generation"=1:1001,"mean"=colMeans(sim[,2:1002]))
heritas50_facevo$lci<-heritas50_facevo$mean-apply(sim[,2:1002],2,sd)/sqrt(max(sim$modelrun))*qt(.975,max(sim$modelrun)-1)
heritas50_facevo$uci<-heritas50_facevo$mean+apply(sim[,2:1002],2,sd)/sqrt(max(sim$modelrun))*qt(.975,max(sim$modelrun)-1)



#Couple Population Analyses#

#Load in the couple populations
setwd(choose.dir(caption="Navigate to folder containing HERIT AS (h = .50) couple populations"))

#Collapse data into a single list
data<-lapply(list.files(),read.csv)

#Tag each model run with a unique modelrun number
for(i in 1:length(data)){
  data[[i]]$modelrun<-i
}

#Collapse all data into a single dataframe
fulldata<-NULL
fulldata<-do.call(rbind.fill,data)
fulldata<-fulldata[,-c(1,25:26)]

#Create a copy of this dataframe to overwrite each loop
procdata<-fulldata

#Create a blank vector for real mate value
procdata$selfmv_main<-NA

#Create a blank vector for mate value based on a random preference point
procdata$selfmv_randpref<-NA

#Create a blank vector for mate value based on scrambled traits
procdata$selfmv_randtrait<-NA

#Do the same for partner mate values
procdata$matemv_main<-NA
procdata$matemv_randpref<-NA
procdata$matemv_randtrait<-NA

#And for d-factor scores
procdata$fa_main<-NA
procdata$fa_randpref<-NA
procdata$fa_randtrait<-NA

#Blank dataframes for storing analysis results
heritas50_vardata<-data.frame("modelrun"=1:max(fulldata$modelrun),"fa_main"=NA,"fa_randpref"=NA,"fa_randtrait"=NA)
heritas50_ccadata<-data.frame("modelrun"=rep(1:max(fulldata$modelrun),each=(20*3)),"factor"=rep(c("fa_main","fa_randpref","fa_randtrait"),each=20),"cca"=NA,"loading"=NA)

for(i in 1:max(fulldata$modelrun)){
  
  loopdata<-fulldata[fulldata$modelrun==i,]
  
  #Split into male and female dataframes
  datam<-subset(loopdata,loopdata$sex==1)
  dataf<-subset(loopdata,loopdata$sex==0)
  
  #Compute the average male and female preferences
  maleprefs<-colMeans(datam[,11:20],na.rm=T)
  femaleprefs<-colMeans(dataf[,11:20],na.rm=T)
  
  #Make a vector of random preferences for each sex
  randmaleprefs<-runif(10,1,7)
  randfemaleprefs<-runif(10,1,7)
  
  #Make dataframes with scrambled traits for each sex
  datamrand<-datam
  datamrand[,1:10]<-apply(datamrand[,1:10],2,function(x) sample(x))
  
  datafrand<-dataf
  datafrand[,1:10]<-apply(datafrand[,1:10],2,function(x) sample(x))
  
  
  
  ###Mate Values###
  #Males#
  #Compute self mate value for males
  datam$selfmv_main<-apply(datam[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_main<-10*(-1*(datam$selfmv_main)+maxeuc)/maxeuc
  
  #Compute self mate value based on a random preference point
  datam$selfmv_randpref<-apply(datam[,1:10],1,function(x) dist(rbind(x,randfemaleprefs)))
  datam$selfmv_randpref<-10*(-1*(datam$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  datam$selfmv_randtrait<-apply(datamrand[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_randtrait<-10*(-1*(datam$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value for males
  datam$matemv_main<-apply(dataf[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_main<-10*(-1*(datam$matemv_main)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  datam$matemv_randpref<-apply(dataf[,1:10],1,function(x) dist(rbind(x,randmaleprefs)))
  datam$matemv_randpref<-10*(-1*(datam$matemv_randpref)+maxeuc)/maxeuc
  
  #Compute partner mate value using scrambled traits
  datam$matemv_randtrait<-apply(datafrand[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_randtrait<-10*(-1*(datam$matemv_randtrait)+maxeuc)/maxeuc
  
  
  #Females#
  #Compute self mate value for females
  dataf$selfmv_main<-apply(dataf[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_main<-10*(-1*(dataf$selfmv_main)+maxeuc)/maxeuc
  
  #Compute self mate value based on a random preference point
  dataf$selfmv_randpref<-apply(dataf[,1:10],1,function(x) dist(rbind(x,randmaleprefs)))
  dataf$selfmv_randpref<-10*(-1*(dataf$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  dataf$selfmv_randtrait<-apply(datafrand[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_randtrait<-10*(-1*(dataf$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value for females
  dataf$matemv_main<-apply(datam[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_main<-10*(-1*(dataf$matemv_main)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  dataf$matemv_randpref<-apply(datam[,1:10],1,function(x) dist(rbind(x,randfemaleprefs)))
  dataf$matemv_randpref<-10*(-1*(dataf$matemv_randpref)+maxeuc)/maxeuc
  
  #Compute partner mate value using scrambled traits
  dataf$matemv_randtrait<-apply(datamrand[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_randtrait<-10*(-1*(dataf$matemv_randtrait)+maxeuc)/maxeuc
  
  
  
  ###Desirabilities###
  #Create dataframes of desirabilities for participants
  selfdesm_main<-data.frame(t(apply(datam[,1:10],1,function(x) 6-abs(x-femaleprefs))))
  selfdesm_randpref<-data.frame(t(apply(datam[,1:10],1,function(x) 6-abs(x-randfemaleprefs))))
  selfdesm_randtrait<-data.frame(t(apply(datamrand[,1:10],1,function(x) 6-abs(x-femaleprefs))))
  
  selfdesf_main<-data.frame(t(apply(dataf[,1:10],1,function(x) 6-abs(x-maleprefs))))
  selfdesf_randpref<-data.frame(t(apply(dataf[,1:10],1,function(x) 6-abs(x-randmaleprefs))))
  selfdesf_randtrait<-data.frame(t(apply(datafrand[,1:10],1,function(x) 6-abs(x-maleprefs))))
  
  #Factor analysis of desirabilities
  mfa_main<-principal(selfdesm_main,factors=1,scores=T)
  mfa_randpref<-principal(selfdesm_randpref,factors=1,scores=T)
  mfa_randtrait<-principal(selfdesm_randtrait,factors=1,scores=T)
  
  ffa_main<-principal(selfdesf_main,factors=1,scores=T)
  ffa_randpref<-principal(selfdesf_randpref,factors=1,scores=T)
  ffa_randtrait<-principal(selfdesf_randtrait,factors=1,scores=T)
  
  ###Variance Explained by Factors###
  #Proportion of variance explained by d-factors for males
  mvar_main<-mfa_main$Vaccounted[2]
  mvar_randpref<-mfa_randpref$Vaccounted[2]
  mvar_randtrait<-mfa_randtrait$Vaccounted[2]
  
  #Proportion of variance explained by d-factors for females
  fvar_main<-ffa_main$Vaccounted[2]
  fvar_randpref<-ffa_randpref$Vaccounted[2]
  fvar_randtrait<-ffa_randtrait$Vaccounted[2]
  
  
  #Average variance explained across males and females
  var_main<-(mvar_main+fvar_main)/2
  var_randpref<-(mvar_randpref+fvar_randpref)/2
  var_randtrait<-(mvar_randtrait+fvar_randtrait)/2
  
  #Save to the appropriate row in the heritas50_vardata dataframe
  heritas50_vardata$fa_main[i]<-var_main
  heritas50_vardata$fa_randpref[i]<-var_randpref
  heritas50_vardata$fa_randtrait[i]<-var_randtrait
  
  
  ###d-Factor Loading as a Function of Cross-Character Assortment###
  
  #Calculate how much each trait generates cross-character assortment for males
  mcca_main<-sapply(1:ncol(selfdesm_main),function(x) mean(cor(selfdesm_main[,x],selfdesf_main[,-x],use="pairwise.complete.obs")))
  mcca_randpref<-sapply(1:ncol(selfdesm_randpref),function(x) mean(cor(selfdesm_randpref[,x],selfdesf_randpref[,-x],use="pairwise.complete.obs")))
  mcca_randtrait<-sapply(1:ncol(selfdesm_randtrait),function(x) mean(cor(selfdesm_randtrait[,x],selfdesf_randtrait[,-x],use="pairwise.complete.obs")))
  
  #Calculate how much each trait generates cross-character assortment for females
  fcca_main<-sapply(1:ncol(selfdesf_main),function(x) mean(cor(selfdesf_main[,x],selfdesm_main[,-x],use="pairwise.complete.obs")))
  fcca_randpref<-sapply(1:ncol(selfdesf_randpref),function(x) mean(cor(selfdesf_randpref[,x],selfdesm_randpref[,-x],use="pairwise.complete.obs")))
  fcca_randtrait<-sapply(1:ncol(selfdesf_randtrait),function(x) mean(cor(selfdesf_randtrait[,x],selfdesm_randtrait[,-x],use="pairwise.complete.obs")))
  
  #Save cross-character assortment values to the cca dataframe
  heritas50_ccadata$cca[heritas50_ccadata$modelrun==i]<-c(fcca_main,mcca_main,fcca_randpref,mcca_randpref,fcca_randtrait,mcca_randtrait)
  
  #Extract the factor loading for each trait
  heritas50_ccadata$loading[heritas50_ccadata$modelrun==i]<-c(as.numeric(ffa_main$loading),as.numeric(mfa_main$loading),as.numeric(ffa_randpref$loading),as.numeric(mfa_randpref$loading),as.numeric(ffa_randtrait$loading),as.numeric(mfa_randtrait$loading))
  
  
  
  ###Save Factor Scores###
  datam$fa_main<-as.numeric(mfa_main$scores)
  datam$fa_randpref<-as.numeric(mfa_randpref$scores)
  datam$fa_randtrait<-as.numeric(mfa_randtrait$scores)
  
  dataf$fa_main<-as.numeric(ffa_main$scores)
  dataf$fa_randpref<-as.numeric(ffa_randpref$scores)
  dataf$fa_randtrait<-as.numeric(ffa_randtrait$scores)
  
  #Standardized mate values
  dataf[,25:30]<-apply(dataf[,25:30],2,function(x) scale(x)*2+5)
  datam[,25:30]<-apply(datam[,25:30],2,function(x) scale(x)*2+5)
  
  #Recombine the data
  loopdata<-rbind(dataf,datam)
  
  #Save the data
  procdata[procdata$modelrun==i,]<-loopdata
  
}

#Rename procdata
heritas50_data<-procdata



###Data Analysis###



##Factor Variance Explained##

#Average variance explained across cultures
heritas50_var<-data.frame("factor"=c("fa_main","fa_randpref","fa_randtrait"))
heritas50_var$heritas50_varexp<-c(mean(heritas50_vardata$fa_main),mean(heritas50_vardata$fa_randpref),mean(heritas50_vardata$fa_randtrait))
heritas50_var$lci<-heritas50_var$heritas50_varexp-c((sd(heritas50_vardata$fa_main)/sqrt(nrow(heritas50_vardata)))*qt(.975,nrow(heritas50_vardata)-1),(sd(heritas50_vardata$fa_randpref)/sqrt(nrow(heritas50_vardata)))*qt(.975,nrow(heritas50_vardata)-1),(sd(heritas50_vardata$fa_randtrait)/sqrt(nrow(heritas50_vardata)))*qt(.975,nrow(heritas50_vardata)-1))
heritas50_var$uci<-heritas50_var$heritas50_varexp+c((sd(heritas50_vardata$fa_main)/sqrt(nrow(heritas50_vardata)))*qt(.975,nrow(heritas50_vardata)-1),(sd(heritas50_vardata$fa_randpref)/sqrt(nrow(heritas50_vardata)))*qt(.975,nrow(heritas50_vardata)-1),(sd(heritas50_vardata$fa_randtrait)/sqrt(nrow(heritas50_vardata)))*qt(.975,nrow(heritas50_vardata)-1))




##Factor Scores and Mate Value##

heritas50_facmv_main<-lmer(fa_main~selfmv_main+(1+selfmv_main|modelrun),data=heritas50_data,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))

#Find the overall line of best fit across model runs
heritas50_facmvovlfit_main<-predict(heritas50_facmv_main)

#Find the coefficients for each model run's regression line
heritas50_facmvcoef_main<-c(summary(heritas50_facmv_main)$coefficients[1,1],summary(heritas50_facmv_main)$coefficients[2,1])

#Standardized Slopes#

heritas50_facmv_mainstd<-stdcoef(heritas50_facmv_main)



##Cross-Character Assortment and d-Factor Loading##

heritas50_cca_main<-lmer(loading~cca+(1+cca|modelrun),data=heritas50_ccadata[heritas50_ccadata$factor=="fa_main",],control = lmerControl(optimizer = "bobyqa"))

#Find the overall line of best fit across model runs
heritas50_ccaovlfit_main<-predict(heritas50_cca_main)

#Find the coefficients for each model run's regression line
heritas50_ccacoef_main<-c(summary(heritas50_cca_main)$coefficients[1,1],summary(heritas50_cca_main)$coefficients[2,1])

#Standardized Slopes#

heritas50_cca_mainstd<-stdcoef(heritas50_cca_main)






##HERIT AS (h=.75)##

#Factor Evolution#

#Load factor size data
setwd(choose.dir(caption="Navigate to folder containing HERIT AS (h=.75) factor variance data"))

#Collapse facdata into a single list
data<-lapply(list.files(),read.csv)

#Tag each model run with a unique modelrun number
for(i in 1:length(data)){
  data[[i]]$modelrun<-i
}

#Collapse all data into a single dataframe
sim<-NULL
sim<-do.call(rbind.fill,data)

heritas75_facevo<-data.frame("model"=rep("Assortative",1001),"generation"=1:1001,"mean"=colMeans(sim[,2:1002]))
heritas75_facevo$lci<-heritas75_facevo$mean-apply(sim[,2:1002],2,sd)/sqrt(max(sim$modelrun))*qt(.975,max(sim$modelrun)-1)
heritas75_facevo$uci<-heritas75_facevo$mean+apply(sim[,2:1002],2,sd)/sqrt(max(sim$modelrun))*qt(.975,max(sim$modelrun)-1)



#Couple Population Analyses#

#Load in the couple populations
setwd(choose.dir(caption="Navigate to folder containing HERIT AS (h=.75) couple populations"))

#Collapse data into a single list
data<-lapply(list.files(),read.csv)

#Tag each model run with a unique modelrun number
for(i in 1:length(data)){
  data[[i]]$modelrun<-i
}

#Collapse all data into a single dataframe
fulldata<-NULL
fulldata<-do.call(rbind.fill,data)
fulldata<-fulldata[,-c(1,25:26)]

#Create a copy of this dataframe to overwrite each loop
procdata<-fulldata

#Create a blank vector for real mate value
procdata$selfmv_main<-NA

#Create a blank vector for mate value based on a random preference point
procdata$selfmv_randpref<-NA

#Create a blank vector for mate value based on scrambled traits
procdata$selfmv_randtrait<-NA

#Do the same for partner mate values
procdata$matemv_main<-NA
procdata$matemv_randpref<-NA
procdata$matemv_randtrait<-NA

#And for d-factor scores
procdata$fa_main<-NA
procdata$fa_randpref<-NA
procdata$fa_randtrait<-NA

#Blank dataframes for storing analysis results
heritas75_vardata<-data.frame("modelrun"=1:max(fulldata$modelrun),"fa_main"=NA,"fa_randpref"=NA,"fa_randtrait"=NA)
heritas75_ccadata<-data.frame("modelrun"=rep(1:max(fulldata$modelrun),each=(20*3)),"factor"=rep(c("fa_main","fa_randpref","fa_randtrait"),each=20),"cca"=NA,"loading"=NA)

for(i in 1:max(fulldata$modelrun)){
  
  loopdata<-fulldata[fulldata$modelrun==i,]
  
  #Split into male and female dataframes
  datam<-subset(loopdata,loopdata$sex==1)
  dataf<-subset(loopdata,loopdata$sex==0)
  
  #Compute the average male and female preferences
  maleprefs<-colMeans(datam[,11:20],na.rm=T)
  femaleprefs<-colMeans(dataf[,11:20],na.rm=T)
  
  #Make a vector of random preferences for each sex
  randmaleprefs<-runif(10,1,7)
  randfemaleprefs<-runif(10,1,7)
  
  #Make dataframes with scrambled traits for each sex
  datamrand<-datam
  datamrand[,1:10]<-apply(datamrand[,1:10],2,function(x) sample(x))
  
  datafrand<-dataf
  datafrand[,1:10]<-apply(datafrand[,1:10],2,function(x) sample(x))
  
  
  
  ###Mate Values###
  #Males#
  #Compute self mate value for males
  datam$selfmv_main<-apply(datam[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_main<-10*(-1*(datam$selfmv_main)+maxeuc)/maxeuc
  
  #Compute self mate value based on a random preference point
  datam$selfmv_randpref<-apply(datam[,1:10],1,function(x) dist(rbind(x,randfemaleprefs)))
  datam$selfmv_randpref<-10*(-1*(datam$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  datam$selfmv_randtrait<-apply(datamrand[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_randtrait<-10*(-1*(datam$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value for males
  datam$matemv_main<-apply(dataf[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_main<-10*(-1*(datam$matemv_main)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  datam$matemv_randpref<-apply(dataf[,1:10],1,function(x) dist(rbind(x,randmaleprefs)))
  datam$matemv_randpref<-10*(-1*(datam$matemv_randpref)+maxeuc)/maxeuc
  
  #Compute partner mate value using scrambled traits
  datam$matemv_randtrait<-apply(datafrand[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_randtrait<-10*(-1*(datam$matemv_randtrait)+maxeuc)/maxeuc
  
  
  #Females#
  #Compute self mate value for females
  dataf$selfmv_main<-apply(dataf[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_main<-10*(-1*(dataf$selfmv_main)+maxeuc)/maxeuc
  
  #Compute self mate value based on a random preference point
  dataf$selfmv_randpref<-apply(dataf[,1:10],1,function(x) dist(rbind(x,randmaleprefs)))
  dataf$selfmv_randpref<-10*(-1*(dataf$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  dataf$selfmv_randtrait<-apply(datafrand[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_randtrait<-10*(-1*(dataf$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value for females
  dataf$matemv_main<-apply(datam[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_main<-10*(-1*(dataf$matemv_main)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  dataf$matemv_randpref<-apply(datam[,1:10],1,function(x) dist(rbind(x,randfemaleprefs)))
  dataf$matemv_randpref<-10*(-1*(dataf$matemv_randpref)+maxeuc)/maxeuc
  
  #Compute partner mate value using scrambled traits
  dataf$matemv_randtrait<-apply(datamrand[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_randtrait<-10*(-1*(dataf$matemv_randtrait)+maxeuc)/maxeuc
  
  
  
  ###Desirabilities###
  #Create dataframes of desirabilities for participants
  selfdesm_main<-data.frame(t(apply(datam[,1:10],1,function(x) 6-abs(x-femaleprefs))))
  selfdesm_randpref<-data.frame(t(apply(datam[,1:10],1,function(x) 6-abs(x-randfemaleprefs))))
  selfdesm_randtrait<-data.frame(t(apply(datamrand[,1:10],1,function(x) 6-abs(x-femaleprefs))))
  
  selfdesf_main<-data.frame(t(apply(dataf[,1:10],1,function(x) 6-abs(x-maleprefs))))
  selfdesf_randpref<-data.frame(t(apply(dataf[,1:10],1,function(x) 6-abs(x-randmaleprefs))))
  selfdesf_randtrait<-data.frame(t(apply(datafrand[,1:10],1,function(x) 6-abs(x-maleprefs))))
  
  #Factor analysis of desirabilities
  mfa_main<-principal(selfdesm_main,factors=1,scores=T)
  mfa_randpref<-principal(selfdesm_randpref,factors=1,scores=T)
  mfa_randtrait<-principal(selfdesm_randtrait,factors=1,scores=T)
  
  ffa_main<-principal(selfdesf_main,factors=1,scores=T)
  ffa_randpref<-principal(selfdesf_randpref,factors=1,scores=T)
  ffa_randtrait<-principal(selfdesf_randtrait,factors=1,scores=T)
  
  ###Variance Explained by Factors###
  #Proportion of variance explained by d-factors for males
  mvar_main<-mfa_main$Vaccounted[2]
  mvar_randpref<-mfa_randpref$Vaccounted[2]
  mvar_randtrait<-mfa_randtrait$Vaccounted[2]
  
  #Proportion of variance explained by d-factors for females
  fvar_main<-ffa_main$Vaccounted[2]
  fvar_randpref<-ffa_randpref$Vaccounted[2]
  fvar_randtrait<-ffa_randtrait$Vaccounted[2]
  
  
  #Average variance explained across males and females
  var_main<-(mvar_main+fvar_main)/2
  var_randpref<-(mvar_randpref+fvar_randpref)/2
  var_randtrait<-(mvar_randtrait+fvar_randtrait)/2
  
  #Save to the appropriate row in the heritas75_vardata dataframe
  heritas75_vardata$fa_main[i]<-var_main
  heritas75_vardata$fa_randpref[i]<-var_randpref
  heritas75_vardata$fa_randtrait[i]<-var_randtrait
  
  
  ###d-Factor Loading as a Function of Cross-Character Assortment###
  
  #Calculate how much each trait generates cross-character assortment for males
  mcca_main<-sapply(1:ncol(selfdesm_main),function(x) mean(cor(selfdesm_main[,x],selfdesf_main[,-x],use="pairwise.complete.obs")))
  mcca_randpref<-sapply(1:ncol(selfdesm_randpref),function(x) mean(cor(selfdesm_randpref[,x],selfdesf_randpref[,-x],use="pairwise.complete.obs")))
  mcca_randtrait<-sapply(1:ncol(selfdesm_randtrait),function(x) mean(cor(selfdesm_randtrait[,x],selfdesf_randtrait[,-x],use="pairwise.complete.obs")))
  
  #Calculate how much each trait generates cross-character assortment for females
  fcca_main<-sapply(1:ncol(selfdesf_main),function(x) mean(cor(selfdesf_main[,x],selfdesm_main[,-x],use="pairwise.complete.obs")))
  fcca_randpref<-sapply(1:ncol(selfdesf_randpref),function(x) mean(cor(selfdesf_randpref[,x],selfdesm_randpref[,-x],use="pairwise.complete.obs")))
  fcca_randtrait<-sapply(1:ncol(selfdesf_randtrait),function(x) mean(cor(selfdesf_randtrait[,x],selfdesm_randtrait[,-x],use="pairwise.complete.obs")))
  
  #Save cross-character assortment values to the cca dataframe
  heritas75_ccadata$cca[heritas75_ccadata$modelrun==i]<-c(fcca_main,mcca_main,fcca_randpref,mcca_randpref,fcca_randtrait,mcca_randtrait)
  
  #Extract the factor loading for each trait
  heritas75_ccadata$loading[heritas75_ccadata$modelrun==i]<-c(as.numeric(ffa_main$loading),as.numeric(mfa_main$loading),as.numeric(ffa_randpref$loading),as.numeric(mfa_randpref$loading),as.numeric(ffa_randtrait$loading),as.numeric(mfa_randtrait$loading))
  
  
  
  ###Save Factor Scores###
  datam$fa_main<-as.numeric(mfa_main$scores)
  datam$fa_randpref<-as.numeric(mfa_randpref$scores)
  datam$fa_randtrait<-as.numeric(mfa_randtrait$scores)
  
  dataf$fa_main<-as.numeric(ffa_main$scores)
  dataf$fa_randpref<-as.numeric(ffa_randpref$scores)
  dataf$fa_randtrait<-as.numeric(ffa_randtrait$scores)
  
  #Standardized mate values
  dataf[,25:30]<-apply(dataf[,25:30],2,function(x) scale(x)*2+5)
  datam[,25:30]<-apply(datam[,25:30],2,function(x) scale(x)*2+5)
  
  #Recombine the data
  loopdata<-rbind(dataf,datam)
  
  #Save the data
  procdata[procdata$modelrun==i,]<-loopdata
  
}

#Rename procdata
heritas75_data<-procdata



###Data Analysis###



##Factor Variance Explained##

#Average variance explained across cultures
heritas75_var<-data.frame("factor"=c("fa_main","fa_randpref","fa_randtrait"))
heritas75_var$heritas75_varexp<-c(mean(heritas75_vardata$fa_main),mean(heritas75_vardata$fa_randpref),mean(heritas75_vardata$fa_randtrait))
heritas75_var$lci<-heritas75_var$heritas75_varexp-c((sd(heritas75_vardata$fa_main)/sqrt(nrow(heritas75_vardata)))*qt(.975,nrow(heritas75_vardata)-1),(sd(heritas75_vardata$fa_randpref)/sqrt(nrow(heritas75_vardata)))*qt(.975,nrow(heritas75_vardata)-1),(sd(heritas75_vardata$fa_randtrait)/sqrt(nrow(heritas75_vardata)))*qt(.975,nrow(heritas75_vardata)-1))
heritas75_var$uci<-heritas75_var$heritas75_varexp+c((sd(heritas75_vardata$fa_main)/sqrt(nrow(heritas75_vardata)))*qt(.975,nrow(heritas75_vardata)-1),(sd(heritas75_vardata$fa_randpref)/sqrt(nrow(heritas75_vardata)))*qt(.975,nrow(heritas75_vardata)-1),(sd(heritas75_vardata$fa_randtrait)/sqrt(nrow(heritas75_vardata)))*qt(.975,nrow(heritas75_vardata)-1))




##Factor Scores and Mate Value##

heritas75_facmv_main<-lmer(fa_main~selfmv_main+(1+selfmv_main|modelrun),data=heritas75_data,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))

#Find the overall line of best fit across model runs
heritas75_facmvovlfit_main<-predict(heritas75_facmv_main)

#Find the coefficients for each model run's regression line
heritas75_facmvcoef_main<-c(summary(heritas75_facmv_main)$coefficients[1,1],summary(heritas75_facmv_main)$coefficients[2,1])

#Standardized Slopes#

heritas75_facmv_mainstd<-stdcoef(heritas75_facmv_main)



##Cross-Character Assortment and d-Factor Loading##

heritas75_cca_main<-lmer(loading~cca+(1+cca|modelrun),data=heritas75_ccadata[heritas75_ccadata$factor=="fa_main",],control=lmerControl(optimizer="Nelder_Mead"))

#Find the overall line of best fit across model runs
heritas75_ccaovlfit_main<-predict(heritas75_cca_main)

#Find the coefficients for each model run's regression line
heritas75_ccacoef_main<-c(summary(heritas75_cca_main)$coefficients[1,1],summary(heritas75_cca_main)$coefficients[2,1])

#Standardized Slopes#

heritas75_cca_mainstd<-stdcoef(heritas75_cca_main)






##HERIT RAND (h=.5)##

#Factor Evolution#

#Load factor size data
setwd(choose.dir(caption="Navigate to folder containing HERIT RAND (h=.5) factor variance data"))

#Collapse facdata into a single list
data<-lapply(list.files(),read.csv)

#Tag each model run with a unique modelrun number
for(i in 1:length(data)){
  data[[i]]$modelrun<-i
}

#Collapse all data into a single dataframe
sim<-NULL
sim<-do.call(rbind.fill,data)

heritrand50_facevo<-data.frame("model"=rep("Assortative",1001),"generation"=1:1001,"mean"=colMeans(sim[,2:1002]))
heritrand50_facevo$lci<-heritrand50_facevo$mean-apply(sim[,2:1002],2,sd)/sqrt(max(sim$modelrun))*qt(.975,max(sim$modelrun)-1)
heritrand50_facevo$uci<-heritrand50_facevo$mean+apply(sim[,2:1002],2,sd)/sqrt(max(sim$modelrun))*qt(.975,max(sim$modelrun)-1)



#Couple Population Analyses#

#Load in the couple populations
setwd(choose.dir(caption="Navigate to folder containing HERIT RAND (h=.5) couple populations"))

#Collapse data into a single list
data<-lapply(list.files(),read.csv)

#Tag each model run with a unique modelrun number
for(i in 1:length(data)){
  data[[i]]$modelrun<-i
}

#Collapse all data into a single dataframe
fulldata<-NULL
fulldata<-do.call(rbind.fill,data)
fulldata<-fulldata[,-c(1,25:26)]

#Create a copy of this dataframe to overwrite each loop
procdata<-fulldata

#Create a blank vector for real mate value
procdata$selfmv_main<-NA

#Create a blank vector for mate value based on a random preference point
procdata$selfmv_randpref<-NA

#Create a blank vector for mate value based on scrambled traits
procdata$selfmv_randtrait<-NA

#Do the same for partner mate values
procdata$matemv_main<-NA
procdata$matemv_randpref<-NA
procdata$matemv_randtrait<-NA

#And for d-factor scores
procdata$fa_main<-NA
procdata$fa_randpref<-NA
procdata$fa_randtrait<-NA

#Blank dataframes for storing analysis results
heritrand50_vardata<-data.frame("modelrun"=1:max(fulldata$modelrun),"fa_main"=NA,"fa_randpref"=NA,"fa_randtrait"=NA)
heritrand50_ccadata<-data.frame("modelrun"=rep(1:max(fulldata$modelrun),each=(20*3)),"factor"=rep(c("fa_main","fa_randpref","fa_randtrait"),each=20),"cca"=NA,"loading"=NA)

for(i in 1:max(fulldata$modelrun)){
  
  loopdata<-fulldata[fulldata$modelrun==i,]
  
  #Split into male and female dataframes
  datam<-subset(loopdata,loopdata$sex==1)
  dataf<-subset(loopdata,loopdata$sex==0)
  
  #Compute the average male and female preferences
  maleprefs<-colMeans(datam[,11:20],na.rm=T)
  femaleprefs<-colMeans(dataf[,11:20],na.rm=T)
  
  #Make a vector of random preferences for each sex
  randmaleprefs<-runif(10,1,7)
  randfemaleprefs<-runif(10,1,7)
  
  #Make dataframes with scrambled traits for each sex
  datamrand<-datam
  datamrand[,1:10]<-apply(datamrand[,1:10],2,function(x) sample(x))
  
  datafrand<-dataf
  datafrand[,1:10]<-apply(datafrand[,1:10],2,function(x) sample(x))
  
  
  
  ###Mate Values###
  #Males#
  #Compute self mate value for males
  datam$selfmv_main<-apply(datam[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_main<-10*(-1*(datam$selfmv_main)+maxeuc)/maxeuc
  
  #Compute self mate value based on a random preference point
  datam$selfmv_randpref<-apply(datam[,1:10],1,function(x) dist(rbind(x,randfemaleprefs)))
  datam$selfmv_randpref<-10*(-1*(datam$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  datam$selfmv_randtrait<-apply(datamrand[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_randtrait<-10*(-1*(datam$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value for males
  datam$matemv_main<-apply(dataf[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_main<-10*(-1*(datam$matemv_main)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  datam$matemv_randpref<-apply(dataf[,1:10],1,function(x) dist(rbind(x,randmaleprefs)))
  datam$matemv_randpref<-10*(-1*(datam$matemv_randpref)+maxeuc)/maxeuc
  
  #Compute partner mate value using scrambled traits
  datam$matemv_randtrait<-apply(datafrand[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_randtrait<-10*(-1*(datam$matemv_randtrait)+maxeuc)/maxeuc
  
  
  #Females#
  #Compute self mate value for females
  dataf$selfmv_main<-apply(dataf[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_main<-10*(-1*(dataf$selfmv_main)+maxeuc)/maxeuc
  
  #Compute self mate value based on a random preference point
  dataf$selfmv_randpref<-apply(dataf[,1:10],1,function(x) dist(rbind(x,randmaleprefs)))
  dataf$selfmv_randpref<-10*(-1*(dataf$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  dataf$selfmv_randtrait<-apply(datafrand[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_randtrait<-10*(-1*(dataf$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value for females
  dataf$matemv_main<-apply(datam[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_main<-10*(-1*(dataf$matemv_main)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  dataf$matemv_randpref<-apply(datam[,1:10],1,function(x) dist(rbind(x,randfemaleprefs)))
  dataf$matemv_randpref<-10*(-1*(dataf$matemv_randpref)+maxeuc)/maxeuc
  
  #Compute partner mate value using scrambled traits
  dataf$matemv_randtrait<-apply(datamrand[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_randtrait<-10*(-1*(dataf$matemv_randtrait)+maxeuc)/maxeuc
  
  
  
  ###Desirabilities###
  #Create dataframes of desirabilities for participants
  selfdesm_main<-data.frame(t(apply(datam[,1:10],1,function(x) 6-abs(x-femaleprefs))))
  selfdesm_randpref<-data.frame(t(apply(datam[,1:10],1,function(x) 6-abs(x-randfemaleprefs))))
  selfdesm_randtrait<-data.frame(t(apply(datamrand[,1:10],1,function(x) 6-abs(x-femaleprefs))))
  
  selfdesf_main<-data.frame(t(apply(dataf[,1:10],1,function(x) 6-abs(x-maleprefs))))
  selfdesf_randpref<-data.frame(t(apply(dataf[,1:10],1,function(x) 6-abs(x-randmaleprefs))))
  selfdesf_randtrait<-data.frame(t(apply(datafrand[,1:10],1,function(x) 6-abs(x-maleprefs))))
  
  #Factor analysis of desirabilities
  mfa_main<-principal(selfdesm_main,factors=1,scores=T)
  mfa_randpref<-principal(selfdesm_randpref,factors=1,scores=T)
  mfa_randtrait<-principal(selfdesm_randtrait,factors=1,scores=T)
  
  ffa_main<-principal(selfdesf_main,factors=1,scores=T)
  ffa_randpref<-principal(selfdesf_randpref,factors=1,scores=T)
  ffa_randtrait<-principal(selfdesf_randtrait,factors=1,scores=T)
  
  ###Variance Explained by Factors###
  #Proportion of variance explained by d-factors for males
  mvar_main<-mfa_main$Vaccounted[2]
  mvar_randpref<-mfa_randpref$Vaccounted[2]
  mvar_randtrait<-mfa_randtrait$Vaccounted[2]
  
  #Proportion of variance explained by d-factors for females
  fvar_main<-ffa_main$Vaccounted[2]
  fvar_randpref<-ffa_randpref$Vaccounted[2]
  fvar_randtrait<-ffa_randtrait$Vaccounted[2]
  
  
  #Average variance explained across males and females
  var_main<-(mvar_main+fvar_main)/2
  var_randpref<-(mvar_randpref+fvar_randpref)/2
  var_randtrait<-(mvar_randtrait+fvar_randtrait)/2
  
  #Save to the appropriate row in the heritrand50_vardata dataframe
  heritrand50_vardata$fa_main[i]<-var_main
  heritrand50_vardata$fa_randpref[i]<-var_randpref
  heritrand50_vardata$fa_randtrait[i]<-var_randtrait
  
  
  ###d-Factor Loading as a Function of Cross-Character Assortment###
  
  #Calculate how much each trait generates cross-character assortment for males
  mcca_main<-sapply(1:ncol(selfdesm_main),function(x) mean(cor(selfdesm_main[,x],selfdesf_main[,-x],use="pairwise.complete.obs")))
  mcca_randpref<-sapply(1:ncol(selfdesm_randpref),function(x) mean(cor(selfdesm_randpref[,x],selfdesf_randpref[,-x],use="pairwise.complete.obs")))
  mcca_randtrait<-sapply(1:ncol(selfdesm_randtrait),function(x) mean(cor(selfdesm_randtrait[,x],selfdesf_randtrait[,-x],use="pairwise.complete.obs")))
  
  #Calculate how much each trait generates cross-character assortment for females
  fcca_main<-sapply(1:ncol(selfdesf_main),function(x) mean(cor(selfdesf_main[,x],selfdesm_main[,-x],use="pairwise.complete.obs")))
  fcca_randpref<-sapply(1:ncol(selfdesf_randpref),function(x) mean(cor(selfdesf_randpref[,x],selfdesm_randpref[,-x],use="pairwise.complete.obs")))
  fcca_randtrait<-sapply(1:ncol(selfdesf_randtrait),function(x) mean(cor(selfdesf_randtrait[,x],selfdesm_randtrait[,-x],use="pairwise.complete.obs")))
  
  #Save cross-character assortment values to the cca dataframe
  heritrand50_ccadata$cca[heritrand50_ccadata$modelrun==i]<-c(fcca_main,mcca_main,fcca_randpref,mcca_randpref,fcca_randtrait,mcca_randtrait)
  
  #Extract the factor loading for each trait
  heritrand50_ccadata$loading[heritrand50_ccadata$modelrun==i]<-c(as.numeric(ffa_main$loading),as.numeric(mfa_main$loading),as.numeric(ffa_randpref$loading),as.numeric(mfa_randpref$loading),as.numeric(ffa_randtrait$loading),as.numeric(mfa_randtrait$loading))
  
  
  
  ###Save Factor Scores###
  datam$fa_main<-as.numeric(mfa_main$scores)
  datam$fa_randpref<-as.numeric(mfa_randpref$scores)
  datam$fa_randtrait<-as.numeric(mfa_randtrait$scores)
  
  dataf$fa_main<-as.numeric(ffa_main$scores)
  dataf$fa_randpref<-as.numeric(ffa_randpref$scores)
  dataf$fa_randtrait<-as.numeric(ffa_randtrait$scores)
  
  #Standardized mate values
  dataf[,25:30]<-apply(dataf[,25:30],2,function(x) scale(x)*2+5)
  datam[,25:30]<-apply(datam[,25:30],2,function(x) scale(x)*2+5)
  
  #Recombine the data
  loopdata<-rbind(dataf,datam)
  
  #Save the data
  procdata[procdata$modelrun==i,]<-loopdata
  
}

#Rename procdata
heritrand50_data<-procdata



###Data Analysis###



##Factor Variance Explained##

#Average variance explained across cultures
heritrand50_var<-data.frame("factor"=c("fa_main","fa_randpref","fa_randtrait"))
heritrand50_var$heritrand50_varexp<-c(mean(heritrand50_vardata$fa_main),mean(heritrand50_vardata$fa_randpref),mean(heritrand50_vardata$fa_randtrait))
heritrand50_var$lci<-heritrand50_var$heritrand50_varexp-c((sd(heritrand50_vardata$fa_main)/sqrt(nrow(heritrand50_vardata)))*qt(.975,nrow(heritrand50_vardata)-1),(sd(heritrand50_vardata$fa_randpref)/sqrt(nrow(heritrand50_vardata)))*qt(.975,nrow(heritrand50_vardata)-1),(sd(heritrand50_vardata$fa_randtrait)/sqrt(nrow(heritrand50_vardata)))*qt(.975,nrow(heritrand50_vardata)-1))
heritrand50_var$uci<-heritrand50_var$heritrand50_varexp+c((sd(heritrand50_vardata$fa_main)/sqrt(nrow(heritrand50_vardata)))*qt(.975,nrow(heritrand50_vardata)-1),(sd(heritrand50_vardata$fa_randpref)/sqrt(nrow(heritrand50_vardata)))*qt(.975,nrow(heritrand50_vardata)-1),(sd(heritrand50_vardata$fa_randtrait)/sqrt(nrow(heritrand50_vardata)))*qt(.975,nrow(heritrand50_vardata)-1))




##Factor Scores and Mate Value##

heritrand50_facmv_main<-lmer(fa_main~selfmv_main+(1+selfmv_main|modelrun),data=heritrand50_data,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))

#Find the overall line of best fit across model runs
heritrand50_facmvovlfit_main<-predict(heritrand50_facmv_main)

#Find the coefficients for each model run's regression line
heritrand50_facmvcoef_main<-c(summary(heritrand50_facmv_main)$coefficients[1,1],summary(heritrand50_facmv_main)$coefficients[2,1])

#Standardized Slopes#

heritrand50_facmv_mainstd<-stdcoef(heritrand50_facmv_main)



##Cross-Character Assortment and d-Factor Loading##

heritrand50_cca_main<-lmer(loading~cca+(1|modelrun),data=heritrand50_ccadata[heritrand50_ccadata$factor=="fa_main",])

#Find the overall line of best fit across model runs
heritrand50_ccaovlfit_main<-predict(heritrand50_cca_main)

#Find the coefficients for each model run's regression line
heritrand50_ccacoef_main<-c(summary(heritrand50_cca_main)$coefficients[1,1],summary(heritrand50_cca_main)$coefficients[2,1])

#Standardized Slopes#

heritrand50_cca_mainstd<-stdcoef(heritrand50_cca_main)







##HERIT RAND (h=.75)##

#Factor Evolution#

#Load factor size data
setwd(choose.dir(caption="Navigate to folder containing HERIT RAND (h=.75) factor variance data"))

#Collapse facdata into a single list
data<-lapply(list.files(),read.csv)

#Tag each model run with a unique modelrun number
for(i in 1:length(data)){
  data[[i]]$modelrun<-i
}

#Collapse all data into a single dataframe
sim<-NULL
sim<-do.call(rbind.fill,data)

heritrand75_facevo<-data.frame("model"=rep("Assortative",1001),"generation"=1:1001,"mean"=colMeans(sim[,2:1002]))
heritrand75_facevo$lci<-heritrand75_facevo$mean-apply(sim[,2:1002],2,sd)/sqrt(max(sim$modelrun))*qt(.975,max(sim$modelrun)-1)
heritrand75_facevo$uci<-heritrand75_facevo$mean+apply(sim[,2:1002],2,sd)/sqrt(max(sim$modelrun))*qt(.975,max(sim$modelrun)-1)



#Couple Population Analyses#

#Load in the couple populations
setwd(choose.dir(caption="Navigate to folder containing HERIT RAND (h=.75) couple populations"))

#Collapse data into a single list
data<-lapply(list.files(),read.csv)

#Tag each model run with a unique modelrun number
for(i in 1:length(data)){
  data[[i]]$modelrun<-i
}

#Collapse all data into a single dataframe
fulldata<-NULL
fulldata<-do.call(rbind.fill,data)
fulldata<-fulldata[,-c(1,25:26)]

#Create a copy of this dataframe to overwrite each loop
procdata<-fulldata

#Create a blank vector for real mate value
procdata$selfmv_main<-NA

#Create a blank vector for mate value based on a random preference point
procdata$selfmv_randpref<-NA

#Create a blank vector for mate value based on scrambled traits
procdata$selfmv_randtrait<-NA

#Do the same for partner mate values
procdata$matemv_main<-NA
procdata$matemv_randpref<-NA
procdata$matemv_randtrait<-NA

#And for d-factor scores
procdata$fa_main<-NA
procdata$fa_randpref<-NA
procdata$fa_randtrait<-NA

#Blank dataframes for storing analysis results
heritrand75_vardata<-data.frame("modelrun"=1:max(fulldata$modelrun),"fa_main"=NA,"fa_randpref"=NA,"fa_randtrait"=NA)
heritrand75_ccadata<-data.frame("modelrun"=rep(1:max(fulldata$modelrun),each=(20*3)),"factor"=rep(c("fa_main","fa_randpref","fa_randtrait"),each=20),"cca"=NA,"loading"=NA)

for(i in 1:max(fulldata$modelrun)){
  
  loopdata<-fulldata[fulldata$modelrun==i,]
  
  #Split into male and female dataframes
  datam<-subset(loopdata,loopdata$sex==1)
  dataf<-subset(loopdata,loopdata$sex==0)
  
  #Compute the average male and female preferences
  maleprefs<-colMeans(datam[,11:20],na.rm=T)
  femaleprefs<-colMeans(dataf[,11:20],na.rm=T)
  
  #Make a vector of random preferences for each sex
  randmaleprefs<-runif(10,1,7)
  randfemaleprefs<-runif(10,1,7)
  
  #Make dataframes with scrambled traits for each sex
  datamrand<-datam
  datamrand[,1:10]<-apply(datamrand[,1:10],2,function(x) sample(x))
  
  datafrand<-dataf
  datafrand[,1:10]<-apply(datafrand[,1:10],2,function(x) sample(x))
  
  
  
  ###Mate Values###
  #Males#
  #Compute self mate value for males
  datam$selfmv_main<-apply(datam[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_main<-10*(-1*(datam$selfmv_main)+maxeuc)/maxeuc
  
  #Compute self mate value based on a random preference point
  datam$selfmv_randpref<-apply(datam[,1:10],1,function(x) dist(rbind(x,randfemaleprefs)))
  datam$selfmv_randpref<-10*(-1*(datam$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  datam$selfmv_randtrait<-apply(datamrand[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_randtrait<-10*(-1*(datam$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value for males
  datam$matemv_main<-apply(dataf[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_main<-10*(-1*(datam$matemv_main)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  datam$matemv_randpref<-apply(dataf[,1:10],1,function(x) dist(rbind(x,randmaleprefs)))
  datam$matemv_randpref<-10*(-1*(datam$matemv_randpref)+maxeuc)/maxeuc
  
  #Compute partner mate value using scrambled traits
  datam$matemv_randtrait<-apply(datafrand[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_randtrait<-10*(-1*(datam$matemv_randtrait)+maxeuc)/maxeuc
  
  
  #Females#
  #Compute self mate value for females
  dataf$selfmv_main<-apply(dataf[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_main<-10*(-1*(dataf$selfmv_main)+maxeuc)/maxeuc
  
  #Compute self mate value based on a random preference point
  dataf$selfmv_randpref<-apply(dataf[,1:10],1,function(x) dist(rbind(x,randmaleprefs)))
  dataf$selfmv_randpref<-10*(-1*(dataf$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  dataf$selfmv_randtrait<-apply(datafrand[,1:10],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_randtrait<-10*(-1*(dataf$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value for females
  dataf$matemv_main<-apply(datam[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_main<-10*(-1*(dataf$matemv_main)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  dataf$matemv_randpref<-apply(datam[,1:10],1,function(x) dist(rbind(x,randfemaleprefs)))
  dataf$matemv_randpref<-10*(-1*(dataf$matemv_randpref)+maxeuc)/maxeuc
  
  #Compute partner mate value using scrambled traits
  dataf$matemv_randtrait<-apply(datamrand[,1:10],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_randtrait<-10*(-1*(dataf$matemv_randtrait)+maxeuc)/maxeuc
  
  
  
  ###Desirabilities###
  #Create dataframes of desirabilities for participants
  selfdesm_main<-data.frame(t(apply(datam[,1:10],1,function(x) 6-abs(x-femaleprefs))))
  selfdesm_randpref<-data.frame(t(apply(datam[,1:10],1,function(x) 6-abs(x-randfemaleprefs))))
  selfdesm_randtrait<-data.frame(t(apply(datamrand[,1:10],1,function(x) 6-abs(x-femaleprefs))))
  
  selfdesf_main<-data.frame(t(apply(dataf[,1:10],1,function(x) 6-abs(x-maleprefs))))
  selfdesf_randpref<-data.frame(t(apply(dataf[,1:10],1,function(x) 6-abs(x-randmaleprefs))))
  selfdesf_randtrait<-data.frame(t(apply(datafrand[,1:10],1,function(x) 6-abs(x-maleprefs))))
  
  #Factor analysis of desirabilities
  mfa_main<-principal(selfdesm_main,factors=1,scores=T)
  mfa_randpref<-principal(selfdesm_randpref,factors=1,scores=T)
  mfa_randtrait<-principal(selfdesm_randtrait,factors=1,scores=T)
  
  ffa_main<-principal(selfdesf_main,factors=1,scores=T)
  ffa_randpref<-principal(selfdesf_randpref,factors=1,scores=T)
  ffa_randtrait<-principal(selfdesf_randtrait,factors=1,scores=T)
  
  ###Variance Explained by Factors###
  #Proportion of variance explained by d-factors for males
  mvar_main<-mfa_main$Vaccounted[2]
  mvar_randpref<-mfa_randpref$Vaccounted[2]
  mvar_randtrait<-mfa_randtrait$Vaccounted[2]
  
  #Proportion of variance explained by d-factors for females
  fvar_main<-ffa_main$Vaccounted[2]
  fvar_randpref<-ffa_randpref$Vaccounted[2]
  fvar_randtrait<-ffa_randtrait$Vaccounted[2]
  
  
  #Average variance explained across males and females
  var_main<-(mvar_main+fvar_main)/2
  var_randpref<-(mvar_randpref+fvar_randpref)/2
  var_randtrait<-(mvar_randtrait+fvar_randtrait)/2
  
  #Save to the appropriate row in the heritrand75_vardata dataframe
  heritrand75_vardata$fa_main[i]<-var_main
  heritrand75_vardata$fa_randpref[i]<-var_randpref
  heritrand75_vardata$fa_randtrait[i]<-var_randtrait
  
  
  ###d-Factor Loading as a Function of Cross-Character Assortment###
  
  #Calculate how much each trait generates cross-character assortment for males
  mcca_main<-sapply(1:ncol(selfdesm_main),function(x) mean(cor(selfdesm_main[,x],selfdesf_main[,-x],use="pairwise.complete.obs")))
  mcca_randpref<-sapply(1:ncol(selfdesm_randpref),function(x) mean(cor(selfdesm_randpref[,x],selfdesf_randpref[,-x],use="pairwise.complete.obs")))
  mcca_randtrait<-sapply(1:ncol(selfdesm_randtrait),function(x) mean(cor(selfdesm_randtrait[,x],selfdesf_randtrait[,-x],use="pairwise.complete.obs")))
  
  #Calculate how much each trait generates cross-character assortment for females
  fcca_main<-sapply(1:ncol(selfdesf_main),function(x) mean(cor(selfdesf_main[,x],selfdesm_main[,-x],use="pairwise.complete.obs")))
  fcca_randpref<-sapply(1:ncol(selfdesf_randpref),function(x) mean(cor(selfdesf_randpref[,x],selfdesm_randpref[,-x],use="pairwise.complete.obs")))
  fcca_randtrait<-sapply(1:ncol(selfdesf_randtrait),function(x) mean(cor(selfdesf_randtrait[,x],selfdesm_randtrait[,-x],use="pairwise.complete.obs")))
  
  #Save cross-character assortment values to the cca dataframe
  heritrand75_ccadata$cca[heritrand75_ccadata$modelrun==i]<-c(fcca_main,mcca_main,fcca_randpref,mcca_randpref,fcca_randtrait,mcca_randtrait)
  
  #Extract the factor loading for each trait
  heritrand75_ccadata$loading[heritrand75_ccadata$modelrun==i]<-c(as.numeric(ffa_main$loading),as.numeric(mfa_main$loading),as.numeric(ffa_randpref$loading),as.numeric(mfa_randpref$loading),as.numeric(ffa_randtrait$loading),as.numeric(mfa_randtrait$loading))
  
  
  
  ###Save Factor Scores###
  datam$fa_main<-as.numeric(mfa_main$scores)
  datam$fa_randpref<-as.numeric(mfa_randpref$scores)
  datam$fa_randtrait<-as.numeric(mfa_randtrait$scores)
  
  dataf$fa_main<-as.numeric(ffa_main$scores)
  dataf$fa_randpref<-as.numeric(ffa_randpref$scores)
  dataf$fa_randtrait<-as.numeric(ffa_randtrait$scores)
  
  #Standardized mate values
  dataf[,25:30]<-apply(dataf[,25:30],2,function(x) scale(x)*2+5)
  datam[,25:30]<-apply(datam[,25:30],2,function(x) scale(x)*2+5)
  
  #Recombine the data
  loopdata<-rbind(dataf,datam)
  
  #Save the data
  procdata[procdata$modelrun==i,]<-loopdata
  
}

#Rename procdata
heritrand75_data<-procdata



###Data Analysis###



##Factor Variance Explained##

#Average variance explained across cultures
heritrand75_var<-data.frame("factor"=c("fa_main","fa_randpref","fa_randtrait"))
heritrand75_var$heritrand75_varexp<-c(mean(heritrand75_vardata$fa_main),mean(heritrand75_vardata$fa_randpref),mean(heritrand75_vardata$fa_randtrait))
heritrand75_var$lci<-heritrand75_var$heritrand75_varexp-c((sd(heritrand75_vardata$fa_main)/sqrt(nrow(heritrand75_vardata)))*qt(.975,nrow(heritrand75_vardata)-1),(sd(heritrand75_vardata$fa_randpref)/sqrt(nrow(heritrand75_vardata)))*qt(.975,nrow(heritrand75_vardata)-1),(sd(heritrand75_vardata$fa_randtrait)/sqrt(nrow(heritrand75_vardata)))*qt(.975,nrow(heritrand75_vardata)-1))
heritrand75_var$uci<-heritrand75_var$heritrand75_varexp+c((sd(heritrand75_vardata$fa_main)/sqrt(nrow(heritrand75_vardata)))*qt(.975,nrow(heritrand75_vardata)-1),(sd(heritrand75_vardata$fa_randpref)/sqrt(nrow(heritrand75_vardata)))*qt(.975,nrow(heritrand75_vardata)-1),(sd(heritrand75_vardata$fa_randtrait)/sqrt(nrow(heritrand75_vardata)))*qt(.975,nrow(heritrand75_vardata)-1))




##Factor Scores and Mate Value##

heritrand75_facmv_main<-lmer(fa_main~selfmv_main+(1+selfmv_main|modelrun),data=heritrand75_data,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))

#Find the overall line of best fit across model runs
heritrand75_facmvovlfit_main<-predict(heritrand75_facmv_main)

#Find the coefficients for each model run's regression line
heritrand75_facmvcoef_main<-c(summary(heritrand75_facmv_main)$coefficients[1,1],summary(heritrand75_facmv_main)$coefficients[2,1])

#Standardized Slopes#

heritrand75_facmv_mainstd<-stdcoef(heritrand75_facmv_main)



##Cross-Character Assortment and d-Factor Loading##

heritrand75_cca_main<-lmer(loading~cca+(1|modelrun),data=heritrand75_ccadata[heritrand75_ccadata$factor=="fa_main",],control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

#Find the overall line of best fit across model runs
heritrand75_ccaovlfit_main<-predict(heritrand75_cca_main)

#Find the coefficients for each model run's regression line
heritrand75_ccacoef_main<-c(summary(heritrand75_cca_main)$coefficients[1,1],summary(heritrand75_cca_main)$coefficients[2,1])

#Standardized Slopes#

heritrand75_cca_mainstd<-stdcoef(heritrand75_cca_main)



###Herit Plots###

#Y-axes
facevoyaxis<-expression(paste("Variance Explaned By ",italic("d"),"-Factor"))
mvplotyaxis<-expression(paste(italic("d"),"-Factor Score"))
ccayaxis<-expression(paste(italic("d"),"-Factor Loading"))


##Herit AS 75 Plots##

##Factor Evolution##

heritas75_facevolimits<-aes(ymin=heritas75_facevo$lci,ymax=heritas75_facevo$uci)
heritas75_facevoplot<-qplot(generation,mean,data=heritas75_facevo,ylab=facevoyaxis,xlab="Generation",geom="line")+theme_classic(base_size=12)+geom_errorbar(heritas75_facevolimits,position="dodge",size=I(.5),alpha=I(.15))+annotate("text",x=min(heritas75_facevo$generation)*1.15,y=max(heritas75_facevo$mean)*1.010,size=8,label="A",fontface="bold")+coord_cartesian(ylim=c(min(heritrand75_facevo$mean)*.99,max(heritas75_facevo$mean)*1.01))


###d-Factor Score as a Function of Mate Value###

#Create plotting dataframes

#Plots
heritas75_facmvplot_main<-qplot(selfmv_main,fa_main,color=as.factor(modelrun),data=heritas75_data,xlab="Agent Mate Value",ylab=mvplotyaxis)+geom_line(aes(y=heritas75_facmvovlfit_main))+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(heritas75_data$selfmv_main)*.95,y=max(heritas75_data$fa_main)*.9,fontface="bold",size=8,label="B")+geom_abline(intercept=heritas75_facmvcoef_main[1],slope=heritas75_facmvcoef_main[2],size=2.5,color="black")



###d-Factor Loading as a Function of Cross-Character Assortment###


heritas75_ccaplot_main<-qplot(cca,loading,color=as.factor(modelrun),data=heritas75_ccadata[heritas75_ccadata$factor=="fa_main",],xlab="Degree of Cross-Character Assortment",ylab=ccayaxis)+geom_line(aes(y=heritas75_ccaovlfit_main),size=1)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(heritas75_ccadata[heritas75_ccadata$factor=="fa_main",]$cca)*.95,y=max(heritas75_ccadata[heritas75_ccadata$factor=="fa_main",]$loading)*.95,fontface="bold",size=8,label="C")+geom_abline(intercept=heritas75_ccacoef_main[1],slope=heritas75_ccacoef_main[2],size=2.5,color="black")





##Herit RAND 75 Plots##

##Factor Evolution##

heritrand75_facevolimits<-aes(ymin=heritrand75_facevo$lci,ymax=heritrand75_facevo$uci)
heritrand75_facevoplot<-qplot(generation,mean,data=heritrand75_facevo,ylab=facevoyaxis,xlab="Generation",geom="line")+theme_classic(base_size=12)+geom_errorbar(heritrand75_facevolimits,position="dodge",size=I(.5),alpha=I(.15))+annotate("text",x=min(heritrand75_facevo$generation)*1.15,y=max(heritas75_facevo$mean)*1.010,size=8,label="D",fontface="bold")+coord_cartesian(ylim=c(min(heritrand75_facevo$mean)*.99,max(heritas75_facevo$mean)*1.01))



###d-Factor Score as a Function of Mate Value###

#Create plotting dataframes

#Plots
heritrand75_facmvplot_main<-qplot(selfmv_main,fa_main,color=as.factor(modelrun),data=heritrand75_data,xlab="Agent Mate Value",ylab=mvplotyaxis)+geom_line(aes(y=heritrand75_facmvovlfit_main))+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(heritrand75_data$selfmv_main)*.95,y=max(heritrand75_data$fa_main)*.9,fontface="bold",size=8,label="E")+geom_abline(intercept=heritrand75_facmvcoef_main[1],slope=heritrand75_facmvcoef_main[2],size=2.5,color="black")



###d-Factor Loading as a Function of Cross-Character Assortment###


heritrand75_ccaplot_main<-qplot(cca,loading,color=as.factor(modelrun),data=heritrand75_ccadata[heritrand75_ccadata$factor=="fa_main",],xlab="Degree of Cross-Character Assortment",ylab=ccayaxis)+geom_line(aes(y=heritrand75_ccaovlfit_main),size=1)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(heritrand75_ccadata[heritrand75_ccadata$factor=="fa_main",]$cca)*.95,y=max(heritrand75_ccadata[heritrand75_ccadata$factor=="fa_main",]$loading)*.95,fontface="bold",size=8,label="F")+geom_abline(intercept=heritrand75_ccacoef_main[1],slope=heritrand75_ccacoef_main[2],size=2.5,color="black")





##Herit AS 50 Plots##

##Factor Evolution##

heritas50_facevolimits<-aes(ymin=heritas50_facevo$lci,ymax=heritas50_facevo$uci)
heritas50_facevoplot<-qplot(generation,mean,data=heritas50_facevo,ylab=facevoyaxis,xlab="Generation",geom="line")+theme_classic(base_size=12)+geom_errorbar(heritas50_facevolimits,position="dodge",size=I(.5),alpha=I(.15))+annotate("text",x=min(heritas50_facevo$generation)*1.15,y=max(heritas50_facevo$mean)*1.010,size=8,label="G",fontface="bold")+coord_cartesian(ylim=c(min(heritrand50_facevo$mean)*.99,max(heritas50_facevo$mean)*1.01))



###d-Factor Score as a Function of Mate Value###

#Create plotting dataframes

#Plots
heritas50_facmvplot_main<-qplot(selfmv_main,fa_main,color=as.factor(modelrun),data=heritas50_data,xlab="Agent Mate Value",ylab=mvplotyaxis)+geom_line(aes(y=heritas50_facmvovlfit_main))+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(heritas50_data$selfmv_main)*.95,y=max(heritas50_data$fa_main)*.9,fontface="bold",size=8,label="H")+geom_abline(intercept=heritas50_facmvcoef_main[1],slope=heritas50_facmvcoef_main[2],size=2.5,color="black")



###d-Factor Loading as a Function of Cross-Character Assortment###


heritas50_ccaplot_main<-qplot(cca,loading,color=as.factor(modelrun),data=heritas50_ccadata[heritas50_ccadata$factor=="fa_main",],xlab="Degree of Cross-Character Assortment",ylab=ccayaxis)+geom_line(aes(y=heritas50_ccaovlfit_main),size=1)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(heritas50_ccadata[heritas50_ccadata$factor=="fa_main",]$cca)*.95,y=max(heritas50_ccadata[heritas50_ccadata$factor=="fa_main",]$loading)*.95,fontface="bold",size=8,label="I")+geom_abline(intercept=heritas50_ccacoef_main[1],slope=heritas50_ccacoef_main[2],size=2.5,color="black")





##Herit RAND 50 Plots##

##Factor Evolution##

heritrand50_facevolimits<-aes(ymin=heritrand50_facevo$lci,ymax=heritrand50_facevo$uci)
heritrand50_facevoplot<-qplot(generation,mean,data=heritrand50_facevo,ylab=facevoyaxis,xlab="Generation",geom="line")+theme_classic(base_size=12)+geom_errorbar(heritrand50_facevolimits,position="dodge",size=I(.5),alpha=I(.15))+annotate("text",x=min(heritrand50_facevo$generation)*1.15,y=max(heritas50_facevo$mean)*1.010,size=8,label="J",fontface="bold")+coord_cartesian(ylim=c(min(heritrand50_facevo$mean)*.99,max(heritas50_facevo$mean)*1.01))



###d-Factor Score as a Function of Mate Value###

#Create plotting dataframes

#Plots
heritrand50_facmvplot_main<-qplot(selfmv_main,fa_main,color=as.factor(modelrun),data=heritrand50_data,xlab="Agent Mate Value",ylab=mvplotyaxis)+geom_line(aes(y=heritrand50_facmvovlfit_main))+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(heritrand50_data$selfmv_main)*.95,y=max(heritrand50_data$fa_main)*.9,fontface="bold",size=8,label="K")+geom_abline(intercept=heritrand50_facmvcoef_main[1],slope=heritrand50_facmvcoef_main[2],size=2.5,color="black")



###d-Factor Loading as a Function of Cross-Character Assortment###


heritrand50_ccaplot_main<-qplot(cca,loading,color=as.factor(modelrun),data=heritrand50_ccadata[heritrand50_ccadata$factor=="fa_main",],xlab="Degree of Cross-Character Assortment",ylab=ccayaxis)+geom_line(aes(y=heritrand50_ccaovlfit_main),size=1)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(heritrand50_ccadata[heritrand50_ccadata$factor=="fa_main",]$cca)*.95,y=max(heritrand50_ccadata[heritrand50_ccadata$factor=="fa_main",]$loading)*.95,fontface="bold",size=8,label="L")+geom_abline(intercept=heritrand50_ccacoef_main[1],slope=heritrand50_ccacoef_main[2],size=2.5,color="black")





###Compile Panel Plots###
herit_plot_main<-grid.arrange(heritas75_facevoplot,heritas50_facevoplot,heritas75_facmvplot_main,heritas75_ccaplot_main,heritas50_facmvplot_main,heritas50_ccaplot_main,heritrand75_facevoplot,heritrand50_facevoplot,heritrand75_facmvplot_main,heritrand75_ccaplot_main,heritrand50_facmvplot_main,heritrand50_ccaplot_main,ncol=4,layout_matrix=cbind(c(1,1,2,2),c(3,4,5,6),c(7,7,8,8),c(9,10,11,12)),widths=c(1,.75,1,.75))



if(sf==1){
  ggsave(paste0(figdir,"\\HERIT Plot Main.tiff"),herit_plot_main,width=16,height=8,unit="in",dpi=300)
  
}

