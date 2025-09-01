######Packages######
library(ggplot2)
library(gridExtra)
library(psych)
library(lme4)
library(lmerTest)
library(reshape2)
library(optimx)



set.seed(06102018)





######Composite Traits######
#Load in the data
data<-read.csv(choose.files(caption="Read in QCS Data"))

#Create a blank dataframe for storing desirability predictive powers and loading
ccadata<-data.frame("cca"=rep(0,40),"loading"=rep(0,40))

#Create a blank dataframe for storing bootstrap factor variance explained results
varloop<-data.frame("loop"=1:100,"var_randpref"=0,"var_randtrait"=0)

#Create a blank dataframe for storing bootstrap mate value results
mvloop<-data.frame("loop"=rep(1:100,each=nrow(data)),"CIN"=data$CIN,"selfmv_randpref"=0,"selfmv_randtrait"=0,"fa_randpref"=0,"fa_randtrait"=0)

#Create a blank dataframe for storing bootstrap predictive power results
ccaloop<-data.frame("loop"=rep(1:100,each=40),"cca_randpref"=0,"cca_randtrait"=0,"loading_randpref"=0,"loading_randtrait"=0)

#Determine the maximum possible Euclidean distance, used for scaling
maxeuc<-sqrt(6^2*20)

###Analyzing the QCS As-is###
#Split into male and female dataframes
datam<-subset(data,data$sex==1)
dataf<-subset(data,data$sex==0)

#Compute the average male and female preferences
maleprefs<-colMeans(datam[,14:33],na.rm=T)
femaleprefs<-colMeans(dataf[,14:33],na.rm=T)


##Mate Values##
#Males#
#Compute self mate value for males
datam$selfmv_main<-apply(datam[,74:93],1,function(x) dist(rbind(x,femaleprefs)))
datam$selfmv_main<-10*(-1*(datam$selfmv_main)+maxeuc)/maxeuc

#Compute partner mate value for males
datam$matemv_main<-apply(dataf[,74:93],1,function(x) dist(rbind(x,maleprefs)))
datam$matemv_main<-10*(-1*(datam$matemv_main)+maxeuc)/maxeuc

#Females#
#Compute self mate value for females
dataf$selfmv_main<-apply(dataf[,74:93],1,function(x) dist(rbind(x,maleprefs)))
dataf$selfmv_main<-10*(-1*(dataf$selfmv_main)+maxeuc)/maxeuc

#Compute partner mate value for females
dataf$matemv_main<-apply(datam[,74:93],1,function(x) dist(rbind(x,femaleprefs)))
dataf$matemv_main<-10*(-1*(dataf$matemv_main)+maxeuc)/maxeuc



##Desirabilities##
#Create dataframes of desirabilities for participants
selfdesm_main<-data.frame(t(apply(datam[,74:93],1,function(x) 6-abs(x-femaleprefs))))

selfdesf_main<-data.frame(t(apply(dataf[,74:93],1,function(x) 6-abs(x-maleprefs))))

#Factor analysis of desirabilities
mfa_main<-principal(selfdesm_main,factors=1,scores=T)

ffa_main<-principal(selfdesf_main,factors=1,scores=T)



##Variance Explained by D-Factor##
#Proportion of variance explained by d-factors for males
mvar_main<-mfa_main$Vaccounted[2]

#Proportion of variance explained by d-factors for males
fvar_main<-ffa_main$Vaccounted[2]

#Average variance explained across males and females
var_main<-(mvar_main+fvar_main)/2



##Factor Loading and Cross-Character Assortment##

#Calculate how much each trait generates cross-character assortment for males
mcca_main<-sapply(1:ncol(selfdesm_main),function(x) mean(cor(selfdesm_main[,x],selfdesf_main[,-x],use="pairwise.complete.obs")))


#Calculate how much each trait generates cross-character assortment for females
fcca_main<-sapply(1:ncol(selfdesf_main),function(x) mean(cor(selfdesf_main[,x],selfdesm_main[,-x],use="pairwise.complete.obs")))

#Save cross-character assortment values to the cca dataframe
ccadata$cca<-c(fcca_main,mcca_main)

#Extract the factor loading for each trait
ccadata$loading<-c(as.numeric(ffa_main$loading),as.numeric(mfa_main$loading))


###Save Factor Scores###
datam$fa_main<-as.numeric(mfa_main$scores)

dataf$fa_main<-as.numeric(ffa_main$scores)




#Rename procdata
data<-rbind(dataf,datam)

#Make couple ID a factor
data$CIN<-as.factor(data$CIN)



######Bootstrapping Random-Point and Scrambled Factor Analyses######

for(l in 1:100){
  ###Random Points and Scrambling###
  
  #Create vectors of random preferences for males and females
  randmaleprefs<-runif(20,1,7)
  randfemaleprefs<-runif(20,1,7)
  
  #Make dataframes with scrambled traits for each sex
  datamrand<-datam
  datamrand[,74:93]<-apply(datamrand[,74:93],2,function(x) sample(x))
  
  datafrand<-dataf
  datafrand[,74:93]<-apply(datafrand[,74:93],2,function(x) sample(x))
  
  
  ###Mate Values###
  
  #Males#
  #Compute self mate value based on a random preference point
  datam$selfmv_randpref<-apply(datam[,74:93],1,function(x) dist(rbind(x,randfemaleprefs)))
  datam$selfmv_randpref<-10*(-1*(datam$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  datam$selfmv_randtrait<-apply(datamrand[,74:93],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_randtrait<-10*(-1*(datam$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  datam$matemv_randpref<-apply(dataf[,74:93],1,function(x) dist(rbind(x,randmaleprefs)))
  datam$matemv_randpref<-10*(-1*(datam$matemv_randpref)+maxeuc)/maxeuc
  
  datam$matemv_randtrait<-apply(datafrand[,74:93],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_randtrait<-10*(-1*(datam$matemv_randtrait)+maxeuc)/maxeuc
  
  #Females#
  #Compute self mate value based on a random preference point
  dataf$selfmv_randpref<-apply(dataf[,74:93],1,function(x) dist(rbind(x,randmaleprefs)))
  dataf$selfmv_randpref<-10*(-1*(dataf$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  dataf$selfmv_randtrait<-apply(datafrand[,74:93],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_randtrait<-10*(-1*(dataf$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  dataf$matemv_randpref<-apply(datam[,74:93],1,function(x) dist(rbind(x,randfemaleprefs)))
  dataf$matemv_randpref<-10*(-1*(dataf$matemv_randpref)+maxeuc)/maxeuc
  
  dataf$matemv_randtrait<-apply(datamrand[,74:93],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_randtrait<-10*(-1*(dataf$matemv_randtrait)+maxeuc)/maxeuc
  
  
  
  ###Desirabilities###
  #Create dataframes of desirabilities for participants
  selfdesm_randpref<-data.frame(t(apply(datam[,74:93],1,function(x) 6-abs(x-randfemaleprefs))))
  selfdesm_randtrait<-data.frame(t(apply(datamrand[,74:93],1,function(x) 6-abs(x-femaleprefs))))
  
  selfdesf_randpref<-data.frame(t(apply(dataf[,74:93],1,function(x) 6-abs(x-randmaleprefs))))
  selfdesf_randtrait<-data.frame(t(apply(datafrand[,74:93],1,function(x) 6-abs(x-maleprefs))))
  
  #Factor analysis of desirabilities
  mfa_randpref<-principal(selfdesm_randpref,factors=1,scores=T)
  mfa_randtrait<-principal(selfdesm_randtrait,factors=1,scores=T)
  
  ffa_randpref<-principal(selfdesf_randpref,factors=1,scores=T)
  ffa_randtrait<-principal(selfdesf_randtrait,factors=1,scores=T)
  
  
  
  ###Variance Explained by Factors###
  #Proportion of variance explained by d-factors for males
  mvar_randpref<-mfa_randpref$Vaccounted[2]
  mvar_randtrait<-mfa_randtrait$Vaccounted[2]
  
  #Proportion of variance explained by d-factors for males
  fvar_randpref<-ffa_randpref$Vaccounted[2]
  fvar_randtrait<-ffa_randtrait$Vaccounted[2]
  
  
  #Average variance explained across males and females
  varloop$var_randpref[l]<-(mvar_randpref+fvar_randpref)/2
  varloop$var_randtrait[l]<-(mvar_randtrait+fvar_randtrait)/2
  
  
  ###Factor Loading and Cross-Character Assortment###
  
  #Calculate how much each trait generates cross-character assortment for males
  mcca_randpref<-sapply(1:ncol(selfdesm_randpref),function(x) mean(cor(selfdesm_randpref[,x],selfdesf_randpref[,-x],use="pairwise.complete.obs")))
  mcca_randtrait<-sapply(1:ncol(selfdesm_randtrait),function(x) mean(cor(selfdesm_randtrait[,x],selfdesf_randtrait[,-x],use="pairwise.complete.obs")))
  
  #Calculate how much each trait generates cross-character assortment for females
  fcca_randpref<-sapply(1:ncol(selfdesf_randpref),function(x) mean(cor(selfdesf_randpref[,x],selfdesm_randpref[,-x],use="pairwise.complete.obs")))
  fcca_randtrait<-sapply(1:ncol(selfdesf_randtrait),function(x) mean(cor(selfdesf_randtrait[,x],selfdesm_randtrait[,-x],use="pairwise.complete.obs")))
  
  #Save cca to the ccaloop dataframe
  ccaloop$cca_randpref[ccaloop$loop==l]<-c(fcca_randpref,mcca_randpref)
  ccaloop$cca_randtrait[ccaloop$loop==l]<-c(fcca_randtrait,mcca_randtrait)
  
  ccaloop$loading_randpref[ccaloop$loop==l]<-c(as.numeric(ffa_randpref$loading),as.numeric(mfa_randpref$loading))
  ccaloop$loading_randtrait[ccaloop$loop==l]<-c(as.numeric(ffa_randtrait$loading),as.numeric(mfa_randtrait$loading))
  
  
  ###Save Factor Scores and Mate Values to mvloop Dataframe###
  mvloop$selfmv_randpref[mvloop$loop==l]<-c(dataf$selfmv_randpref,datam$selfmv_randpref)
  mvloop$selfmv_randtrait[mvloop$loop==l]<-c(dataf$selfmv_randtrait,datam$selfmv_randtrait)
  
  mvloop$fa_randpref[mvloop$loop==l]<-c(as.numeric(ffa_randpref$scores),as.numeric(mfa_randpref$scores))
  mvloop$fa_randtrait[mvloop$loop==l]<-c(as.numeric(ffa_randtrait$scores),as.numeric(mfa_randtrait$scores))
  
}






###Data Analysis###

##Factor Variance Explained##

#Average variance explained across cultures
var_comp<-data.frame("factor"=c("fa_main","fa_randpref","fa_randtrait"))
var_comp$varexp<-c(var_main,mean(varloop$var_randpref),mean(varloop$var_randtrait))
var_comp$lci<-var_comp$varexp-c(NA,(sd(varloop$var_randpref)/sqrt(nrow(varloop)))*qt(.975,nrow(varloop)-1),(sd(varloop$var_randtrait)/sqrt(nrow(varloop)))*qt(.975,nrow(varloop)-1))
var_comp$uci<-var_comp$varexp+c(NA,(sd(varloop$var_randpref)/sqrt(nrow(varloop)))*qt(.975,nrow(varloop)-1),(sd(varloop$var_randtrait)/sqrt(nrow(varloop)))*qt(.975,nrow(varloop)-1))


##Factor Scores and Mate Value##
facmv_main_comp<-lmer(fa_main~selfmv_main+(1|CIN),data=data)
facmv_randpref_comp<-lmer(fa_randpref~selfmv_randpref+(1|CIN)+(1+selfmv_randpref|loop),data=mvloop,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))
facmv_randtrait_comp<-lmer(fa_randtrait~selfmv_randtrait+(1|CIN)+(1+selfmv_randtrait|loop),data=mvloop,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

facmv_mainstd_comp<-stdcoef(facmv_main_comp)
facmv_randprefstd_comp<-stdcoef(facmv_randpref_comp)
facmv_randtraitstd_comp<-stdcoef(facmv_randtrait_comp)

facmvfit_randpref_comp<-predict(facmv_randpref_comp,mvloop,allow.new.levels=T)
facmvfit_randtrait_comp<-predict(facmv_randtrait_comp,mvloop,allow.new.levels=T)

mvslopes_randpref_comp<-coef(facmv_randpref_comp)$loop
mvslopes_randtrait_comp<-coef(facmv_randtrait_comp)$loop

facmvcoef_randpref_comp<-c(summary(facmv_randpref_comp)$coefficients[1,1],summary(facmv_randpref_comp)$coefficients[2,1])
facmvcoef_randtrait_comp<-c(summary(facmv_randtrait_comp)$coefficients[1,1],summary(facmv_randtrait_comp)$coefficients[2,1])



##Cross-Character Assortment and Factor Loading
cca_main_comp<-lm(loading~cca,data=ccadata)
cca_randpref_comp<-lmer(loading_randpref~cca_randpref+(1+cca_randpref|loop),data=ccaloop,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))
cca_randtrait_comp<-lmer(loading_randtrait~cca_randtrait+(1+cca_randtrait|loop),data=ccaloop,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

cca_mainstd_comp<-lm(scale(loading)~scale(cca),data=ccadata)
cca_randprefstd_comp<-stdcoef(cca_randpref_comp)
cca_randtraitstd_comp<-stdcoef(cca_randtrait_comp)

loadfit_randpref_comp<-predict(cca_randpref_comp,ccaloop,allow.new.levels=T)
loadfit_randtrait_comp<-predict(cca_randtrait_comp,ccaloop,allow.new.levels=T)

loadslopes_randpref_comp<-coef(cca_randpref_comp)$loop
loadslopes_randtrait_comp<-coef(cca_randtrait_comp)$loop

ccacoef_randpref_comp<-c(summary(cca_randpref_comp)$coefficients[1,1],summary(cca_randpref_comp)$coefficients[2,1])
ccacoef_randtrait_comp<-c(summary(cca_randtrait_comp)$coefficients[1,1],summary(cca_randtrait_comp)$coefficients[2,1])



###Plots###

###d-Factor Score as a Function of Mate Value###

mvplotyaxis<-expression(paste(italic("d"),"-Factor Score"))

#Plots
facmvplot_main_comp<-qplot(selfmv_main,fa_main,data=data,xlab="Participant Mate Value",ylab=mvplotyaxis)+geom_smooth(method="lm",se=F)+theme_classic(base_size=12)+annotate("text",x=4,y=1.25,fontface="bold",size=8,label="A")


###d-Factor Loading as a Function of Predictive Power###
#Plots

predyaxis<-expression(paste(italic("d"),"-Factor loading"))

ccaplot_main_comp<-qplot(cca,loading,data=ccadata,xlab="Degree of Cross-Character Assortment",ylab=predyaxis)+geom_smooth(method="lm",se=F)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=.05,y=.65,fontface="bold",size=8,label="B")










######Mate Ratings Traits######

#Load in the data
data<-read.csv(choose.files("Load in the QCS Data"))

# data<-data[complete.cases(data[,54:73]),]

#Create a blank dataframe for storing desirability predictive powers and loading
ccadata<-data.frame("cca"=rep(0,40),"loading"=rep(0,40))

#Create a blank dataframe for storing bootstrap factor variance explained results
varloop<-data.frame("loop"=1:100,"var_randpref"=0,"var_randtrait"=0)

#Create a blank dataframe for storing bootstrap mate value results
mvloop<-data.frame("loop"=rep(1:100,each=nrow(data)),"CIN"=data$CIN,"selfmv_randpref"=0,"selfmv_randtrait"=0,"fa_randpref"=0,"fa_randtrait"=0)

#Create a blank dataframe for storing bootstrap predictive power results
ccaloop<-data.frame("loop"=rep(1:100,each=40),"cca_randpref"=0,"cca_randtrait"=0,"loading_randpref"=0,"loading_randtrait"=0)

#Determine the maximum possible Euclidean distance, used for scaling
maxeuc<-sqrt(6^2*20)

###Analyzing the QCS As-is###
#Split into male and female dataframes
datam<-subset(data,data$sex==1)
dataf<-subset(data,data$sex==0)

#Compute the average male and female preferences
maleprefs<-colMeans(datam[,14:33],na.rm=T)
femaleprefs<-colMeans(dataf[,14:33],na.rm=T)


##Mate Values##
#Males#
#Compute self mate value for males
datam$selfmv_main<-apply(datam[,54:73],1,function(x) dist(rbind(x,femaleprefs)))
datam$selfmv_main<-10*(-1*(datam$selfmv_main)+maxeuc)/maxeuc

#Compute partner mate value for males
datam$matemv_main<-apply(dataf[,54:73],1,function(x) dist(rbind(x,maleprefs)))
datam$matemv_main<-10*(-1*(datam$matemv_main)+maxeuc)/maxeuc

#Females#
#Compute self mate value for females
dataf$selfmv_main<-apply(dataf[,54:73],1,function(x) dist(rbind(x,maleprefs)))
dataf$selfmv_main<-10*(-1*(dataf$selfmv_main)+maxeuc)/maxeuc

#Compute partner mate value for females
dataf$matemv_main<-apply(datam[,54:73],1,function(x) dist(rbind(x,femaleprefs)))
dataf$matemv_main<-10*(-1*(dataf$matemv_main)+maxeuc)/maxeuc



##Desirabilities##
#Create dataframes of desirabilities for participants
selfdesm_main<-data.frame(t(apply(datam[,54:73],1,function(x) 6-abs(x-femaleprefs))))

selfdesf_main<-data.frame(t(apply(dataf[,54:73],1,function(x) 6-abs(x-maleprefs))))

#Factor analysis of desirabilities
mfa_main<-principal(selfdesm_main,factors=1,scores=T)

ffa_main<-principal(selfdesf_main,factors=1,scores=T)



##Variance Explained by D-Factor##
#Proportion of variance explained by d-factors for males
mvar_main<-mfa_main$Vaccounted[2]

#Proportion of variance explained by d-factors for males
fvar_main<-ffa_main$Vaccounted[2]

#Average variance explained across males and females
var_main<-(mvar_main+fvar_main)/2



##Factor Loading as a Function of Predictive Power##

#Calculate how predictive each desirability is of partner mate value for males
mcca_main<-apply(selfdesm_main,2,function(x) cor(x,datam$matemv_main,use="pairwise.complete.obs"))

#Calculate how predictive each desirability is of partner mate value for males
fcca_main<-apply(selfdesf_main,2,function(x) cor(x,dataf$matemv_main,use="pairwise.complete.obs"))

#Save cca to the pred dataframe
ccadata$cca<-c(fcca_main,mcca_main)

#Extract the factor loading for each trait
ccadata$loading<-c(as.numeric(ffa_main$loading),as.numeric(mfa_main$loading))

ccadata$intercor<-c(diag(cor(dataf[,34:53],dataf[,54:73],use="pairwise.complete.obs")),diag(cor(datam[,34:53],dataf[,54:73],use="pairwise.complete.obs")))


###Save Factor Scores###
datam$fa_main<-as.numeric(mfa_main$scores)

dataf$fa_main<-as.numeric(ffa_main$scores)




#Rename procdata
data<-rbind(dataf,datam)

#Make couple ID a factor
data$CIN<-as.factor(data$CIN)



######Bootstrapping Random-Point and Scrambled Factor Analyses######

for(l in 1:100){
  ###Random Points and Scrambling###
  
  #Create vectors of random preferences for males and females
  randmaleprefs<-runif(20,1,7)
  randfemaleprefs<-runif(20,1,7)
  
  #Make dataframes with scrambled traits for each sex
  datamrand<-datam
  datamrand[,54:73]<-apply(datamrand[,54:73],2,function(x) sample(x))
  
  datafrand<-dataf
  datafrand[,54:73]<-apply(datafrand[,54:73],2,function(x) sample(x))
  
  
  ###Mate Values###
  
  #Males#
  #Compute self mate value based on a random preference point
  datam$selfmv_randpref<-apply(datam[,54:73],1,function(x) dist(rbind(x,randfemaleprefs)))
  datam$selfmv_randpref<-10*(-1*(datam$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  datam$selfmv_randtrait<-apply(datamrand[,54:73],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_randtrait<-10*(-1*(datam$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  datam$matemv_randpref<-apply(dataf[,54:73],1,function(x) dist(rbind(x,randmaleprefs)))
  datam$matemv_randpref<-10*(-1*(datam$matemv_randpref)+maxeuc)/maxeuc
  
  datam$matemv_randtrait<-apply(datafrand[,54:73],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_randtrait<-10*(-1*(datam$matemv_randtrait)+maxeuc)/maxeuc
  
  #Females#
  #Compute self mate value based on a random preference point
  dataf$selfmv_randpref<-apply(dataf[,54:73],1,function(x) dist(rbind(x,randmaleprefs)))
  dataf$selfmv_randpref<-10*(-1*(dataf$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  dataf$selfmv_randtrait<-apply(datafrand[,54:73],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_randtrait<-10*(-1*(dataf$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  dataf$matemv_randpref<-apply(datam[,54:73],1,function(x) dist(rbind(x,randfemaleprefs)))
  dataf$matemv_randpref<-10*(-1*(dataf$matemv_randpref)+maxeuc)/maxeuc
  
  dataf$matemv_randtrait<-apply(datamrand[,54:73],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_randtrait<-10*(-1*(dataf$matemv_randtrait)+maxeuc)/maxeuc
  
  
  
  ###Desirabilities###
  #Create dataframes of desirabilities for participants
  selfdesm_randpref<-data.frame(t(apply(datam[,54:73],1,function(x) 6-abs(x-randfemaleprefs))))
  selfdesm_randtrait<-data.frame(t(apply(datamrand[,54:73],1,function(x) 6-abs(x-femaleprefs))))
  
  selfdesf_randpref<-data.frame(t(apply(dataf[,54:73],1,function(x) 6-abs(x-randmaleprefs))))
  selfdesf_randtrait<-data.frame(t(apply(datafrand[,54:73],1,function(x) 6-abs(x-maleprefs))))
  
  #Factor analysis of desirabilities
  mfa_randpref<-principal(selfdesm_randpref,factors=1,scores=T)
  mfa_randtrait<-principal(selfdesm_randtrait,factors=1,scores=T)
  
  ffa_randpref<-principal(selfdesf_randpref,factors=1,scores=T)
  ffa_randtrait<-principal(selfdesf_randtrait,factors=1,scores=T)
  
  
  
  ###Variance Explained by Factors###
  #Proportion of variance explained by d-factors for males
  mvar_randpref<-mfa_randpref$Vaccounted[2]
  mvar_randtrait<-mfa_randtrait$Vaccounted[2]
  
  #Proportion of variance explained by d-factors for males
  fvar_randpref<-ffa_randpref$Vaccounted[2]
  fvar_randtrait<-ffa_randtrait$Vaccounted[2]
  
  
  #Average variance explained across males and females
  varloop$var_randpref[l]<-(mvar_randpref+fvar_randpref)/2
  varloop$var_randtrait[l]<-(mvar_randtrait+fvar_randtrait)/2
  
  
  ###Factor Loading as a Function of Predictive Power###
  
  #Calculate how predictive each desirability is of partner mate value for males
  mcca_randpref<-apply(selfdesm_randpref,2,function(x) cor(x,datam$matemv_randpref,use="pairwise.complete.obs"))
  mcca_randtrait<-apply(selfdesm_randtrait,2,function(x) cor(x,datam$matemv_randtrait,use="pairwise.complete.obs"))
  
  #Calculate how predictive each desirability is of partner mate value for males
  fcca_randpref<-apply(selfdesf_randpref,2,function(x) cor(x,dataf$matemv_randpref,use="pairwise.complete.obs"))
  fcca_randtrait<-apply(selfdesf_randtrait,2,function(x) cor(x,dataf$matemv_randtrait,use="pairwise.complete.obs"))
  
  #Save cca to the ccaloop dataframe
  ccaloop$cca_randpref[ccaloop$loop==l]<-c(fcca_randpref,mcca_randpref)
  ccaloop$cca_randtrait[ccaloop$loop==l]<-c(fcca_randtrait,mcca_randtrait)
  
  ccaloop$loading_randpref[ccaloop$loop==l]<-c(as.numeric(ffa_randpref$loading),as.numeric(mfa_randpref$loading))
  ccaloop$loading_randtrait[ccaloop$loop==l]<-c(as.numeric(ffa_randtrait$loading),as.numeric(mfa_randtrait$loading))
  
  
  ###Save Factor Scores and Mate Values to mvloop Dataframe###
  mvloop$selfmv_randpref[mvloop$loop==l]<-c(dataf$selfmv_randpref,datam$selfmv_randpref)
  mvloop$selfmv_randtrait[mvloop$loop==l]<-c(dataf$selfmv_randtrait,datam$selfmv_randtrait)
  
  mvloop$fa_randpref[mvloop$loop==l]<-c(as.numeric(ffa_randpref$scores),as.numeric(mfa_randpref$scores))
  mvloop$fa_randtrait[mvloop$loop==l]<-c(as.numeric(ffa_randtrait$scores),as.numeric(mfa_randtrait$scores))
  
}






###Data Analysis###

##Factor Variance Explained##

#Average variance explained across cultures
var_mr<-data.frame("factor"=c("fa_main","fa_randpref","fa_randtrait"))
var_mr$varexp<-c(var_main,mean(varloop$var_randpref),mean(varloop$var_randtrait))
var_mr$lci<-var_mr$varexp-c(NA,(sd(varloop$var_randpref)/sqrt(nrow(varloop)))*qt(.975,nrow(varloop)-1),(sd(varloop$var_randtrait)/sqrt(nrow(varloop)))*qt(.975,nrow(varloop)-1))
var_mr$uci<-var_mr$varexp+c(NA,(sd(varloop$var_randpref)/sqrt(nrow(varloop)))*qt(.975,nrow(varloop)-1),(sd(varloop$var_randtrait)/sqrt(nrow(varloop)))*qt(.975,nrow(varloop)-1))


##Factor Scores and Mate Value##
facmv_main_mr<-lmer(fa_main~selfmv_main+(1|CIN),data=data)

facmv_mainstd_mr<-lmer(scale(fa_main)~scale(selfmv_main)+(1|CIN),data=data)



##Cross Character Assortment and Factor Loading##
cca_main_mr<-lm(loading~cca,data=ccadata)

cca_mainstd_mr<-lm(scale(loading)~scale(cca),data=ccadata)


##Inter-Rater Agreement and Factor Loading##
cca_inter<-lm(intercor~loading,data=ccadata)



###Plots###

###d-Factor Score as a Function of Mate Value###

mvplotyaxis<-expression(paste(italic("d"),"-Factor Score"))

#Plots
facmvplot_main_mr<-qplot(selfmv_main,fa_main,data=data,xlab="Participant Mate Value",ylab=mvplotyaxis)+geom_smooth(method="lm",se=F)+theme_classic(base_size=12)+annotate("text",x=4,y=1.1,fontface="bold",size=8,label="C")


###d-Factor Loading as a Function of Predictive Power###
#Plots

predyaxis<-expression(paste(italic("d"),"-Factor loading"))

ccaplot_main_mr<-qplot(cca,loading,data=ccadata,xlab="Degree of Cross-Character Assortment",ylab=predyaxis)+geom_smooth(method="lm",se=F)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=-.1,y=.60,fontface="bold",size=8,label="D")


plot_main<-grid.arrange(facmvplot_main_comp,facmvplot_main_mr,ccaplot_main_comp,ccaplot_main_mr,ncol=2,layout_matrix=cbind(c(1,2),c(3,4)))

interplot_main<-qplot(loading,intercor,data=ccadata,xlab="Desirability Factor Loading",ylab="Inter-Rater Agreement Correlation")+geom_smooth(method="lm",se=F)+theme_classic(base_size=20)
