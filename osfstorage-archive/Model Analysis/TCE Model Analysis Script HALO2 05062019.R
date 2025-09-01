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




###HALO2 OTHER###

#HALO2 analyses based exclusively on other reports

#Couple Population Analyses#

#Load in the couple populations
setwd(choose.dir(caption="Navigate to folder containing HALO2 couple population data"))

#Collapse data into a single list
data<-lapply(list.files(),read.csv)

#Tag each model run with a unique modelrun number
for(i in 1:length(data)){
  data[[i]]$modelrun<-i
}

#Collapse all data into a single dataframe
fulldata<-NULL
fulldata<-do.call(rbind.fill,data)
fulldata<-fulldata[,-1]

#Create a copy of this dataframe to overwrite each loop
procdata<-fulldata

#Create a blank vector for real mate value
procdata$selfmv_main<-NA


#Do the same for partner mate values
procdata$matemv_main<-NA

#And for d-factor scores
procdata$fa_main<-NA

#Blank dataframes for storing analysis results
halo2other_vardata<-data.frame("modelrun"=1:max(fulldata$modelrun),"fa_main"=NA)
halo2other_ccadata<-data.frame("modelrun"=rep(1:max(fulldata$modelrun),each=20),"cca"=NA,"loading"=NA,"intercor"=NA)

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
  
  ###Mate Values###
  #Males#
  #Compute self mate value for males (based on female ratings)
  datam$selfmv_main<-apply(dataf[,21:30],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_main<-10*(-1*(datam$selfmv_main)+maxeuc)/maxeuc
  
  #Compute partner mate value for males (based on male "self-report")
  datam$matemv_main<-apply(datam[,21:30],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_main<-10*(-1*(datam$matemv_main)+maxeuc)/maxeuc
  
  
  #Females#
  #Compute self mate value for females
  dataf$selfmv_main<-apply(datam[,21:30],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_main<-10*(-1*(dataf$selfmv_main)+maxeuc)/maxeuc
  
  #Compute partner mate value for females (based on female "self-report")
  dataf$matemv_main<-apply(dataf[,21:30],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_main<-10*(-1*(dataf$matemv_main)+maxeuc)/maxeuc
  
  ###Desirabilities###
  #Create dataframes of desirabilities for participants
  selfdesm_main<-data.frame(t(apply(dataf[,21:30],1,function(x) 6-abs(x-femaleprefs))))
  selfdesf_main<-data.frame(t(apply(datam[,21:30],1,function(x) 6-abs(x-maleprefs))))
  
  matedesm_main<-data.frame(t(apply(datam[,21:30],1,function(x) 6-abs(x-maleprefs))))
  matedesf_main<-data.frame(t(apply(dataf[,21:30],1,function(x) 6-abs(x-femaleprefs))))
  
  #Factor analysis of desirabilities
  mfa_main<-principal(selfdesm_main,factors=1,scores=T)
  ffa_main<-principal(selfdesf_main,factors=1,scores=T)
  
  ###Variance Explained by Factors###
  #Proportion of variance explained by d-factors for males
  mvar_main<-mfa_main$Vaccounted[2]
  
  #Proportion of variance explained by d-factors for females
  fvar_main<-ffa_main$Vaccounted[2]
  
  
  #Average variance explained across males and females
  var_main<-(mvar_main+fvar_main)/2
  
  #Save to the appropriate row in the halo2other_vardata dataframe
  halo2other_vardata$fa_main[i]<-var_main
  
  
  ###d-Factor Loading as a Function of Cross-Character Assortment###
  
  #Calculate how much each trait generates cross-character assortment for males
  mcca_main<-sapply(1:ncol(selfdesm_main),function(x) mean(cor(selfdesm_main[,x],matedesm_main[,-x],use="pairwise.complete.obs")))
  
  #Calculate how much each trait generates cross-character assortment for females
  fcca_main<-sapply(1:ncol(selfdesf_main),function(x) mean(cor(selfdesf_main[,x],matedesf_main[,-x],use="pairwise.complete.obs")))
  
  #Calculate the inter-rater agreement correlation for all traits
  mintercor<-diag(cor(datam[,1:10],dataf[,21:30]))
  fintercor<-diag(cor(dataf[,1:10],datam[,21:30]))
  
  #Save cross-character assortment values to the cca dataframe
  halo2other_ccadata$cca[halo2other_ccadata$modelrun==i]<-c(fcca_main,mcca_main)
  halo2other_ccadata$intercor[halo2other_ccadata$modelrun==i]<-c(fintercor,mintercor)
  
  
  #Extract the factor loading for each trait
  halo2other_ccadata$loading[halo2other_ccadata$modelrun==i]<-c(as.numeric(ffa_main$loading),as.numeric(mfa_main$loading))
  
  
  
  ###Save Factor Scores###
  datam$fa_main<-as.numeric(mfa_main$scores)
  
  dataf$fa_main<-as.numeric(ffa_main$scores)
  
  #Standardized mate values
  dataf[,35:36]<-apply(dataf[,35:36],2,function(x) scale(x)*2+5)
  datam[,35:36]<-apply(datam[,35:36],2,function(x) scale(x)*2+5)
  
  #Recombine the data
  loopdata<-rbind(dataf,datam)
  
  #Save the data
  procdata[procdata$modelrun==i,]<-loopdata
  
}

#Rename procdata
halo2other_data<-procdata



###Data Analysis###


##Assortative Mating for MV##

#Correlation between partner mate values across model runs
mvcors<-sapply(unique(halo2other_data$modelrun),function(x) cor(halo2other_data[halo2other_data$modelrun==x & halo2other_data$sex==0,]$selfmv_main,halo2other_data[halo2other_data$modelrun==x & halo2other_data$sex==0,]$matemv_main))

#Calculate mean and 95% CIs on assortative mating correlations
mvcordata<-data.frame("mean"=mean(mvcors))
mvcordata$lci<-mvcordata$mean-(sd(mvcors)/sqrt(max(halo2other_data$modelrun))*qt(.975,max(halo2other_data$modelrun)-1))
mvcordata$uci<-mvcordata$mean+(sd(mvcors)/sqrt(max(halo2other_data$modelrun))*qt(.975,max(halo2other_data$modelrun)-1))


##Factor Variance Explained##

#Average variance explained across model runs
halo2other_var<-data.frame("factor"="fa_main")
halo2other_var$halo2other_varexp<-mean(halo2other_vardata$fa_main)
halo2other_var$lci<-halo2other_var$halo2other_varexp-(sd(halo2other_vardata$fa_main)/sqrt(nrow(halo2other_vardata)))*qt(.975,nrow(halo2other_vardata)-1)
halo2other_var$uci<-halo2other_var$halo2other_varexp+(sd(halo2other_vardata$fa_main)/sqrt(nrow(halo2other_vardata)))*qt(.975,nrow(halo2other_vardata)-1)




##Factor Scores and Mate Value##

halo2other_facmv_main<-lmer(fa_main~selfmv_main+(1+selfmv_main|modelrun),data=halo2other_data,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

#Find the overall line of best fit across model runs
halo2other_facmvovlfit_main<-predict(halo2other_facmv_main)

#Find the coefficients for each model run's regression line
halo2other_facmvcoef_main<-c(summary(halo2other_facmv_main)$coefficients[1,1],summary(halo2other_facmv_main)$coefficients[2,1])

#Standardized Slopes#

halo2other_facmv_mainstd<-stdcoef(halo2other_facmv_main)



##Cross-Character Assortment and d-Factor Loading##

halo2other_cca_main<-lmer(loading~cca+(1+cca|modelrun),data=halo2other_ccadata,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))

#Find the overall line of best fit across model runs
halo2other_ccaovlfit_main<-predict(halo2other_cca_main)

#Find the coefficients for each model run's regression line
halo2other_ccacoef_main<-c(summary(halo2other_cca_main)$coefficients[1,1],summary(halo2other_cca_main)$coefficients[2,1])

#Standardized Slopes#

halo2other_cca_mainstd<-stdcoef(halo2other_cca_main)



##d-Factor Loading and Inter-Rater Agreement##

halo2other_inter_main<-lmer(intercor~loading+(1+loading|modelrun),data=halo2other_ccadata,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))

#Find the overall line of best fit across model runs
halo2other_interfit_main<-predict(halo2other_inter_main)

#Find the coefficients for each model run's regression line
halo2other_intercoef_main<-c(summary(halo2other_inter_main)$coefficients[1,1],summary(halo2other_inter_main)$coefficients[2,1])



###Plots###

##Variance Explained by d-Factor##

#Create the plotting dataframe
halo2other_varplotdata<-melt(halo2other_vardata,id.vars=c("modelrun"))

varyaxis<-expression(paste("Variance Explaned By",italic(" d"),"-Factor"))

#Plots
halo2other_varplot_main<-qplot(as.factor(modelrun),value,data=halo2other_varplotdata[halo2other_varplotdata$variable=="fa_main",],fill=as.factor(modelrun),xlab="Modelrun",ylab=varyaxis,geom="blank")+geom_bar(stat="identity",position="dodge")+theme_classic(base_size=12)+theme(legend.position="none")+scale_fill_discrete(name="Factor")+geom_hline(yintercept=halo2other_var[1,2],size=2)+annotate("text",x=5,y=max(halo2other_varplotdata[halo2other_varplotdata$variable=="fa_main",]$value)*1.05,fontface="bold",size=8,label="A")



###d-Factor Score as a Function of Mate Value###

#Create plotting dataframes
mvplotyaxis<-expression(paste(italic("d"),"-Factor Score"))

#Plots
halo2other_facmvplot_main<-qplot(selfmv_main,fa_main,color=as.factor(modelrun),data=halo2other_data,xlab="Agent Mate Value",ylab=mvplotyaxis)+geom_line(aes(y=halo2other_facmvovlfit_main))+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(halo2other_data$selfmv_main)*.95,y=max(halo2other_data$fa_main)*.9,fontface="bold",size=8,label="B")+geom_abline(intercept=halo2other_facmvcoef_main[1],slope=halo2other_facmvcoef_main[2],size=2.5,color="black")



###d-Factor Loading as a Function of Cross-Character Assortment###

ccayaxis<-expression(paste(italic("d"),"-Factor Loading"))

halo2other_ccaplot_main<-qplot(cca,loading,color=as.factor(modelrun),data=halo2other_ccadata,xlab="Degree of Cross-Character Assortment",ylab=ccayaxis)+geom_line(aes(y=halo2other_ccaovlfit_main),size=1)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(halo2other_ccadata$cca)*.95,y=max(halo2other_ccadata$loading)*.95,fontface="bold",size=8,label="C")+geom_abline(intercept=halo2other_ccacoef_main[1],slope=halo2other_ccacoef_main[2],size=2.5,color="black")


###Compile Panel Plots###
halo2other_plot_main<-grid.arrange(halo2other_varplot_main,halo2other_facmvplot_main,halo2other_ccaplot_main,ncol=2,layout_matrix=cbind(c(1,1),c(2,3)),widths=c(1,.75))
halo2other_loadinghist<-qplot(halo2other_ccadata$loading,xlab="Desirability Factor Loading",ylab="Frequency",binwidth=.05,fill=I("coral"),color=I("black"))+theme_classic(base_size=20)
halo2other_interplot_main<-qplot(loading,intercor,color=as.factor(modelrun),data=halo2other_ccadata,xlab="Desirability Factor Loading",ylab="Inter-Rater Agreement Correlation")+geom_line(aes(y=halo2other_interfit_main),size=1)+theme_classic(base_size=20)+theme(legend.position="none")+geom_abline(intercept=halo2other_intercoef_main[1],slope=halo2other_intercoef_main[2],size=2.5,color="black")



if(sf==1){
  ggsave(paste0(figdir,"\\Halo2 Plot Main.tiff"),halo2other_plot_main,width=10,height=5,unit="in",dpi=300)
  ggsave(paste0(figdir,"\\Halo2 Plot Loading Hist.tiff"),halo2other_loadinghist,width=10,height=5,unit="in",dpi=300)
  ggsave(paste0(figdir,"\\Halo2 Plot Inter.tiff"),halo2other_interplot_main,width=10,height=5,unit="in",dpi=300)
  
}


