######Packages######
library(ggplot2)
library(gridExtra)
library(psych)
library(lme4)
library(lmerTest)
library(reshape2)
library(optimx)



######Functions######
stdcoef <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}



######Seed######
#Set a seed to make the random-point and scrambled factor analysis results reproducible
set.seed(06102018)



######Data######
#Load in the data
data<-choose.files(caption="Load in the Newlywed Data")

#Remove two participants for whom we don't have trait ratings
data<-data[-c(90,197),]

#Create a blank dataframe for storing cross-character assortment and loading
ccadata<-data.frame("cca"=rep(0,80),"loading"=rep(0,80))

#Create a blank dataframe for storing bootstrap factor variance explained results
varloop<-data.frame("loop"=1:100,"var_randpref"=0,"var_randtrait"=0)

#Create a blank dataframe for storing bootstrap mate value results
mvloop<-data.frame("loop"=rep(1:100,each=nrow(data)),"CIN"=data$CIN,"selfmv_randpref"=0,"selfmv_randtrait"=0,"fa_randpref"=0,"fa_randtrait"=0)

#Create a blank dataframe for storing bootstrap predictive power results
ccaloop<-data.frame("loop"=rep(1:100,each=80),"cca_randpref"=0,"cca_randtrait"=0,"loading_randpref"=0,"loading_randtrait"=0)

#Determine the maximum possible Euclidean distance, used for scaling
maxeuc<-sqrt(6^2*40)

###Analyzing the Newlywed Data As-is###
#Split into male and female dataframes
datam<-subset(data,data$sex==1)
dataf<-subset(data,data$sex==0)

#Compute the average male and female preferences
maleprefs<-colMeans(datam[,4:43],na.rm=T)
femaleprefs<-colMeans(dataf[,4:43],na.rm=T)



##Mate Values##
#Males#
#Compute self mate value for males
datam$selfmv_main<-apply(datam[,44:83],1,function(x) dist(rbind(x,femaleprefs)))
datam$selfmv_main<-10*(-1*(datam$selfmv_main)+maxeuc)/maxeuc

#Compute partner mate value for males
datam$matemv_main<-apply(dataf[,44:83],1,function(x) dist(rbind(x,maleprefs)))
datam$matemv_main<-10*(-1*(datam$matemv_main)+maxeuc)/maxeuc

#Females#
#Compute self mate value for females
dataf$selfmv_main<-apply(dataf[,44:83],1,function(x) dist(rbind(x,maleprefs)))
dataf$selfmv_main<-10*(-1*(dataf$selfmv_main)+maxeuc)/maxeuc

#Compute partner mate value for females
dataf$matemv_main<-apply(datam[,44:83],1,function(x) dist(rbind(x,femaleprefs)))
dataf$matemv_main<-10*(-1*(dataf$matemv_main)+maxeuc)/maxeuc



##Desirabilities##
#Create dataframes of desirabilities for participants
selfdesm_main<-data.frame(t(apply(datam[,44:83],1,function(x) 6-abs(x-femaleprefs))))

selfdesf_main<-data.frame(t(apply(dataf[,44:83],1,function(x) 6-abs(x-maleprefs))))

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
  randmaleprefs<-runif(40,1,7)
  randfemaleprefs<-runif(40,1,7)
  
  #Make dataframes with scrambled traits for each sex
  datamrand<-datam
  datamrand[,44:83]<-apply(datamrand[,44:83],2,function(x) sample(x))
  
  datafrand<-dataf
  datafrand[,44:83]<-apply(datafrand[,44:83],2,function(x) sample(x))
  
  
  ###Mate Values###
  
  #Males#
  #Compute self mate value based on a random preference point
  datam$selfmv_randpref<-apply(datam[,44:83],1,function(x) dist(rbind(x,randfemaleprefs)))
  datam$selfmv_randpref<-10*(-1*(datam$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  datam$selfmv_randtrait<-apply(datamrand[,44:83],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_randtrait<-10*(-1*(datam$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  datam$matemv_randpref<-apply(dataf[,44:83],1,function(x) dist(rbind(x,randmaleprefs)))
  datam$matemv_randpref<-10*(-1*(datam$matemv_randpref)+maxeuc)/maxeuc
  
  datam$matemv_randtrait<-apply(datafrand[,44:83],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_randtrait<-10*(-1*(datam$matemv_randtrait)+maxeuc)/maxeuc
  
  #Females#
  #Compute self mate value based on a random preference point
  dataf$selfmv_randpref<-apply(dataf[,44:83],1,function(x) dist(rbind(x,randmaleprefs)))
  dataf$selfmv_randpref<-10*(-1*(dataf$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  dataf$selfmv_randtrait<-apply(datafrand[,44:83],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_randtrait<-10*(-1*(dataf$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  dataf$matemv_randpref<-apply(datam[,44:83],1,function(x) dist(rbind(x,randfemaleprefs)))
  dataf$matemv_randpref<-10*(-1*(dataf$matemv_randpref)+maxeuc)/maxeuc
  
  dataf$matemv_randtrait<-apply(datamrand[,44:83],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_randtrait<-10*(-1*(dataf$matemv_randtrait)+maxeuc)/maxeuc
  
  
  
  ###Desirabilities###
  #Create dataframes of desirabilities for participants
  selfdesm_randpref<-data.frame(t(apply(datam[,44:83],1,function(x) 6-abs(x-randfemaleprefs))))
  selfdesm_randtrait<-data.frame(t(apply(datamrand[,44:83],1,function(x) 6-abs(x-femaleprefs))))
  
  selfdesf_randpref<-data.frame(t(apply(dataf[,44:83],1,function(x) 6-abs(x-randmaleprefs))))
  selfdesf_randtrait<-data.frame(t(apply(datafrand[,44:83],1,function(x) 6-abs(x-maleprefs))))
  
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
var<-data.frame("factor"=c("fa_main","fa_randpref","fa_randtrait"))
var$varexp<-c(var_main,mean(varloop$var_randpref),mean(varloop$var_randtrait))
var$lci<-var$varexp-c(NA,(sd(varloop$var_randpref)/sqrt(nrow(varloop)))*qt(.975,nrow(varloop)-1),(sd(varloop$var_randtrait)/sqrt(nrow(varloop)))*qt(.975,nrow(varloop)-1))
var$uci<-var$varexp+c(NA,(sd(varloop$var_randpref)/sqrt(nrow(varloop)))*qt(.975,nrow(varloop)-1),(sd(varloop$var_randtrait)/sqrt(nrow(varloop)))*qt(.975,nrow(varloop)-1))


##Factor Scores and Mate Value##
facmv_main<-lmer(fa_main~selfmv_main+(1|CIN),data=data)
facmv_randpref<-lmer(fa_randpref~selfmv_randpref+(1|CIN)+(1+selfmv_randpref|loop),data=mvloop)
facmv_randtrait<-lmer(fa_randtrait~selfmv_randtrait+(1|CIN)+(1+selfmv_randtrait|loop),data=mvloop)

facmv_mainstd<-stdcoef(facmv_main)
facmv_randprefstd<-stdcoef(facmv_randpref)
facmv_randtraitstd<-stdcoef(facmv_randtrait)

facmvfit_randpref<-predict(facmv_randpref)
facmvfit_randtrait<-predict(facmv_randtrait)

mvslopes_randpref<-coef(facmv_randpref)$loop
mvslopes_randtrait<-coef(facmv_randtrait)$loop

facmvcoef_randpref<-c(summary(facmv_randpref)$coefficients[1,1],summary(facmv_randpref)$coefficients[2,1])
facmvcoef_randtrait<-c(summary(facmv_randtrait)$coefficients[1,1],summary(facmv_randtrait)$coefficients[2,1])



##Factor Loading and Cross-Character Assortment
cca_main<-lm(loading~cca,data=ccadata)
cca_randpref<-lmer(loading_randpref~cca_randpref+(1+cca_randpref|loop),data=ccaloop,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))
cca_randtrait<-lmer(loading_randtrait~cca_randtrait+(1+cca_randtrait|loop),data=ccaloop)

cca_mainstd<-lm(scale(loading)~scale(cca),data=ccadata)
cca_randprefstd<-stdcoef(cca_randpref)
cca_randtraitstd<-stdcoef(cca_randtrait)

loadfit_randpref<-predict(cca_randpref)
loadfit_randtrait<-predict(cca_randtrait)

loadslopes_randpref<-coef(cca_randpref)$loop
loadslopes_randtrait<-coef(cca_randtrait)$loop

ccacoef_randpref<-c(summary(cca_randpref)$coefficients[1,1],summary(cca_randpref)$coefficients[2,1])
ccacoef_randtrait<-c(summary(cca_randtrait)$coefficients[1,1],summary(cca_randtrait)$coefficients[2,1])



###Plots###

###d-Factor Score as a Function of Mate Value###

#Plots
facmvplot_main<-qplot(selfmv_main,fa_main,data=data,xlab="Participant Mate Value",ylab=mvplotyaxis)+geom_smooth(method="lm",se=F)+theme_classic(base_size=12)+annotate("text",x=6.5,y=1.6,fontface="bold",size=8,label="A")

facmvplot_randpref<-qplot(selfmv_randpref,fa_randpref,color=as.factor(loop),data=mvloop,xlab="Participant Mate Value",ylab=mvplotyaxis)+geom_line(aes(y=facmvfit_randpref))+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=4.6,y=3.5,fontface="bold",size=8,label="A")+geom_abline(intercept=facmvcoef_randpref[1],slope=facmvcoef_randpref[2],size=2.5,color="black")
facmvplot_randtrait<-qplot(selfmv_randtrait,fa_randtrait,color=as.factor(loop),data=mvloop,xlab="Participant Mate Value",ylab=mvplotyaxis)+geom_line(aes(y=facmvfit_randtrait))+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=7,y=3.5,fontface="bold",size=8,label="A")+geom_abline(intercept=facmvcoef_randtrait[1],slope=facmvcoef_randtrait[2],size=2.5,color="black")

###d-Factor Loading as a Function of Predictive Power###
#Plots

predyaxis<-expression(paste(italic("d"),"-Factor loading"))

cca_main<-qplot(cca,loading,data=ccadata,xlab="Degree of Cross-Character Assortment",ylab=predyaxis)+geom_smooth(method="lm",se=F)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=-.05,y=.95,fontface="bold",size=8,label="B")
cca_randpref<-qplot(cca_randpref,loading_randpref,color=as.factor(loop),data=ccaloop,xlab="Degree of Cross-Character Assortment",ylab=predyaxis)+geom_line(aes(y=loadfit_randpref),size=1)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=-.20,y=.95,fontface="bold",size=8,label="B")+geom_abline(intercept=ccacoef_randpref[1],slope=ccacoef_randpref[2],size=2.5,color="black")
cca_randtrait<-qplot(cca_randtrait,loading_randtrait,color=as.factor(loop),data=ccaloop,xlab="Degree of Cross-Character Assortment",ylab=predyaxis)+geom_line(aes(y=loadfit_randtrait),size=1)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=-.05,y=.75,fontface="bold",size=8,label="B")+geom_abline(intercept=ccacoef_randtrait[1],slope=ccacoef_randtrait[2],size=2.5,color="black")


plot_main<-grid.arrange(facmvplot_main,cca_main,ncol=2)
plot_randpref<-grid.arrange(facmvplot_randpref,cca_randpref,ncol=2)
plot_randtrait<-grid.arrange(facmvplot_randtrait,cca_randtrait,ncol=2)
