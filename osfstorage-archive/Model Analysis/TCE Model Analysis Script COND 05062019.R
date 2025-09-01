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






###COND###


#Factor Evolution#

#Load factor size data
setwd(choose.dir(caption="Navigate to folder containing COND factor variance data"))

#Collapse facdata into a single list
data<-lapply(list.files(),read.csv)

#Tag each model run with a unique modelrun number
for(i in 1:length(data)){
  data[[i]]$modelrun<-i
}

#Collapse all data into a single dataframe
sim<-NULL
sim<-do.call(rbind.fill,data)

cond_facevo<-data.frame("model"=rep("Assortative",1001),"generation"=1:1001,"mean"=colMeans(sim[,2:1002]))
cond_facevo$lci<-cond_facevo$mean-apply(sim[,2:1002],2,sd)/sqrt(max(sim$modelrun))*qt(.975,max(sim$modelrun)-1)
cond_facevo$uci<-cond_facevo$mean+apply(sim[,2:1002],2,sd)/sqrt(max(sim$modelrun))*qt(.975,max(sim$modelrun)-1)



#Couple Population Analyses#

#Load in the couple populations
setwd(choose.dir(caption="Navigate to folder containing COND couple populations"))

#Collapse data into a single list
data<-lapply(list.files(),read.csv)

#Tag each model run with a unique modelrun number
for(i in 1:length(data)){
  data[[i]]$modelrun<-i
}

#Collapse all data into a single dataframe
fulldata<-NULL
fulldata<-do.call(rbind.fill,data)
fulldata<-fulldata[,-c(1,26:27)]

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
cond_vardata<-data.frame("modelrun"=1:max(fulldata$modelrun),"fa_main"=NA,"fa_randpref"=NA,"fa_randtrait"=NA)
cond_ccadata<-data.frame("modelrun"=rep(1:max(fulldata$modelrun),each=(20*3)),"factor"=rep(c("fa_main","fa_randpref","fa_randtrait"),each=20),"cca"=NA,"loading"=NA)

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
  
  #Save to the appropriate row in the cond_vardata dataframe
  cond_vardata$fa_main[i]<-var_main
  cond_vardata$fa_randpref[i]<-var_randpref
  cond_vardata$fa_randtrait[i]<-var_randtrait
  
  
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
  cond_ccadata$cca[cond_ccadata$modelrun==i]<-c(fcca_main,mcca_main,fcca_randpref,mcca_randpref,fcca_randtrait,mcca_randtrait)
  
  #Extract the factor loading for each trait
  cond_ccadata$loading[cond_ccadata$modelrun==i]<-c(as.numeric(ffa_main$loading),as.numeric(mfa_main$loading),as.numeric(ffa_randpref$loading),as.numeric(mfa_randpref$loading),as.numeric(ffa_randtrait$loading),as.numeric(mfa_randtrait$loading))
  
  
  
  ###Save Factor Scores###
  datam$fa_main<-as.numeric(mfa_main$scores)
  datam$fa_randpref<-as.numeric(mfa_randpref$scores)
  datam$fa_randtrait<-as.numeric(mfa_randtrait$scores)
  
  dataf$fa_main<-as.numeric(ffa_main$scores)
  dataf$fa_randpref<-as.numeric(ffa_randpref$scores)
  dataf$fa_randtrait<-as.numeric(ffa_randtrait$scores)
  
  #Standardized mate values
  dataf[,26:31]<-apply(dataf[,26:31],2,function(x) scale(x)*2+5)
  datam[,26:31]<-apply(datam[,26:31],2,function(x) scale(x)*2+5)
  
  #Recombine the data
  loopdata<-rbind(dataf,datam)
  
  #Save the data
  procdata[procdata$modelrun==i,]<-loopdata
  
}

#Rename procdata
cond_data<-procdata



###Data Analysis###



##Factor Variance Explained##

#Average variance explained across cultures
cond_var<-data.frame("factor"=c("fa_main","fa_randpref","fa_randtrait"))
cond_var$cond_varexp<-c(mean(cond_vardata$fa_main),mean(cond_vardata$fa_randpref),mean(cond_vardata$fa_randtrait))
cond_var$lci<-cond_var$cond_varexp-c((sd(cond_vardata$fa_main)/sqrt(nrow(cond_vardata)))*qt(.975,nrow(cond_vardata)-1),(sd(cond_vardata$fa_randpref)/sqrt(nrow(cond_vardata)))*qt(.975,nrow(cond_vardata)-1),(sd(cond_vardata$fa_randtrait)/sqrt(nrow(cond_vardata)))*qt(.975,nrow(cond_vardata)-1))
cond_var$uci<-cond_var$cond_varexp+c((sd(cond_vardata$fa_main)/sqrt(nrow(cond_vardata)))*qt(.975,nrow(cond_vardata)-1),(sd(cond_vardata$fa_randpref)/sqrt(nrow(cond_vardata)))*qt(.975,nrow(cond_vardata)-1),(sd(cond_vardata$fa_randtrait)/sqrt(nrow(cond_vardata)))*qt(.975,nrow(cond_vardata)-1))




##Factor Scores and Mate Value##

cond_facmv_main<-lmer(fa_main~selfmv_main+(1+selfmv_main|modelrun),data=cond_data,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

#Find the overall line of best fit across model runs
cond_facmvovlfit_main<-predict(cond_facmv_main)

#Find the coefficients for each model run's regression line
cond_facmvcoef_main<-c(summary(cond_facmv_main)$coefficients[1,1],summary(cond_facmv_main)$coefficients[2,1])

#Standardized Slopes#

cond_facmv_mainstd<-stdcoef(cond_facmv_main)



##Cross-Character Assortment and d-Factor Loading##

cond_cca_main<-lmer(loading~cca+(1+cca|modelrun),data=cond_ccadata[cond_ccadata$factor=="fa_main",],control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))

#Find the overall line of best fit across model runs
cond_ccaovlfit_main<-predict(cond_cca_main)

#Find the coefficients for each model run's regression line
cond_ccacoef_main<-c(summary(cond_cca_main)$coefficients[1,1],summary(cond_cca_main)$coefficients[2,1])

#Standardized Slopes#

cond_cca_mainstd<-stdcoef(cond_cca_main)


###Plots###

##Factor Evolution##

yaxis<-expression(paste("Variance Explaned By ",italic("d"),"-Factor"))

limits<-aes(ymin=cond_facevo$lci,ymax=cond_facevo$uci)
cond_facevoplot<-qplot(generation,mean,data=cond_facevo,ylab=yaxis,xlab="Generation",geom="line")+theme_classic(base_size=12)+geom_errorbar(limits,position="dodge",size=I(.5),alpha=I(.15))+annotate("text",x=min(cond_facevo$generation)*1.15,y=max(cond_facevo$mean)*1.010,size=8,label="A",fontface="bold")+coord_cartesian(ylim=c(min(cond_facevo$mean)*.99,max(cond_facevo$mean)*1.01))



##Variance Explained by d-Factor##

#Create the plotting dataframe
cond_varplotdata<-melt(cond_vardata,id.vars=c("modelrun"))

varyaxis<-expression(paste("Variance Explaned By",italic(" d"),"-Factor"))

#Plots
cond_varplot_main<-qplot(as.factor(modelrun),value,data=cond_varplotdata[cond_varplotdata$variable=="fa_main",],fill=as.factor(modelrun),xlab="Modelrun",ylab=varyaxis,geom="blank")+geom_bar(stat="identity",position="dodge")+theme_classic(base_size=12)+theme(legend.position="none")+scale_fill_discrete(name="Factor")+geom_hline(yintercept=cond_var[1,2],size=2)+annotate("text",x=5,y=max(cond_varplotdata[cond_varplotdata$variable=="fa_main",]$value)*1.05,fontface="bold",size=8,label="A")



###d-Factor Score as a Function of Mate Value###

#Create plotting dataframes
mvplotyaxis<-expression(paste(italic("d"),"-Factor Score"))

#Plots
cond_facmvplot_main<-qplot(selfmv_main,fa_main,color=as.factor(modelrun),data=cond_data,xlab="Agent Mate Value",ylab=mvplotyaxis)+geom_line(aes(y=cond_facmvovlfit_main))+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(cond_data$selfmv_main)*.95,y=max(cond_data$fa_main)*.9,fontface="bold",size=8,label="B")+geom_abline(intercept=cond_facmvcoef_main[1],slope=cond_facmvcoef_main[2],size=2.5,color="black")



###d-Factor Loading as a Function of Cross-Character Assortment###

ccayaxis<-expression(paste(italic("d"),"-Factor Loading"))

cond_ccaplot_main<-qplot(cca,loading,color=as.factor(modelrun),data=cond_ccadata[cond_ccadata$factor=="fa_main",],xlab="Degree of Cross-Character Assortment",ylab=ccayaxis)+geom_line(aes(y=cond_ccaovlfit_main),size=1)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(cond_ccadata[cond_ccadata$factor=="fa_main",]$cca)*.95,y=max(cond_ccadata[cond_ccadata$factor=="fa_main",]$loading)*.95,fontface="bold",size=8,label="C")+geom_abline(intercept=cond_ccacoef_main[1],slope=cond_ccacoef_main[2],size=2.5,color="black")



###Compile Panel Plots###
cond_plot_main<-grid.arrange(cond_facevoplot,cond_facmvplot_main,cond_ccaplot_main,ncol=2,layout_matrix=cbind(c(1,1),c(2,3)),widths=c(1,.75))




if(sf==1){
  ggsave(paste0(figdir,"\\Cond Plot Main.tiff"),cond_plot_main,width=10,height=5,unit="in",dpi=300)

}

