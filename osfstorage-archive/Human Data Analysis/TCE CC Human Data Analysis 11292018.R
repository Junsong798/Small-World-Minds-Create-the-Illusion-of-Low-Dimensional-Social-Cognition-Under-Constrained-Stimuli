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
set.seed(06082018)



######Data######
data<-read.csv(choose.files(caption="Load in Human Cross-Cultural Data"))

#Proportion of participants from community populations
pcommunity<-sum(data$sample==2,na.rm=T)/sum(data$sample==2 | data$sample==1,na.rm=T)

#Age information
agedata<-summary(data$age)

#Proportion of people in a commited romantic relationship
prelstat<-sum(data$relstat==1,na.rm=T)/nrow(data)

#Frequency and proportion of each relationship type
allrels<-data.frame("reltype"=c("Single","Dating","Engaged","Married","Divorced","Widowed"),"freq"=as.numeric(table(data$relstat2)),"prop"=as.numeric(table(data$relstat2))/sum(as.numeric(table(data$relstat2))))

#Frequency and proportion of just commited relationship types
crels<-data.frame("reltype"=c("Dating","Engaged","Married"),"freq"=as.numeric(table(data$relstat2))[2:4],"prop"=as.numeric(table(data$relstat2))[2:4]/sum(as.numeric(table(data$relstat2))[2:4]))

#Countrywise sample sizes
csample<-data.frame(table(data$country))
colnames(csample)<-c("country","n")


#A separate dataframe to be looped through for calculations
procdata<-data

#Create a blank vector for real mate value
procdata$selfmv_main<-NA

#Create a blank vector for mate value based on a random preference point
procdata$selfmv_randpref<-NA

#Create a blank vector for mate value based on scrambled traits
procdata$selfmv_randtrait<-NA


procdata$matemv_main<-NA
procdata$matemv_randpref<-NA
procdata$matemv_randtrait<-NA

procdata$fa_main<-NA
procdata$fa_randpref<-NA
procdata$fa_randtrait<-NA

procdata$des_intelligence<-NA
procdata$des_kindness<-NA
procdata$des_health<-NA
procdata$des_physatt<-NA
procdata$des_resources<-NA

#Blank dataframes for storing analysis results
vardata<-data.frame("CIN"=unique(data$CIN),"fa_main"=NA,"fa_randpref"=NA,"fa_randtrait"=NA)
ccadata<-data.frame("CIN"=rep(unique(data$CIN),each=(10*3)),"factor"=rep(c("fa_main","fa_randpref","fa_randtrait"),each=10),"cca"=NA,"loading"=NA)


maxeuc<-sqrt(6^2*5)


for(i in 1:max(data$CIN)){
  ###Processing###
  #Select just the data from the ith country
  loopdata<-subset(data,data$CIN==i)
  
  #Split into male and female dataframes
  datam<-subset(loopdata,loopdata$sex==1)
  dataf<-subset(loopdata,loopdata$sex==0)
  
  #Compute the average male and female preferences
  maleprefs<-colMeans(datam[,15:19],na.rm=T)
  femaleprefs<-colMeans(dataf[,15:19],na.rm=T)
  
  #Make a vector of random preferences for each sex
  randmaleprefs<-runif(5,1,7)
  randfemaleprefs<-runif(5,1,7)
  
  #Make dataframes with scrambled traits for each sex
  datamrand<-datam
  datamrand[,20:29]<-apply(datamrand[,20:29],2,function(x) sample(x))
  
  datafrand<-dataf
  datafrand[,20:29]<-apply(datafrand[,20:29],2,function(x) sample(x))
  
  
  
  ###Mate Values###
  #Males#
  #Compute self mate value for males
  datam$selfmv_main<-apply(datam[,20:24],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_main<-10*(-1*(datam$selfmv_main)+maxeuc)/maxeuc
  
  #Compute self mate value based on a random preference point
  datam$selfmv_randpref<-apply(datam[,20:24],1,function(x) dist(rbind(x,randfemaleprefs)))
  datam$selfmv_randpref<-10*(-1*(datam$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  datam$selfmv_randtrait<-apply(datamrand[,20:24],1,function(x) dist(rbind(x,femaleprefs)))
  datam$selfmv_randtrait<-10*(-1*(datam$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value for males
  datam$matemv_main<-apply(datam[,25:29],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_main<-10*(-1*(datam$matemv_main)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  datam$matemv_randpref<-apply(datam[,25:29],1,function(x) dist(rbind(x,randmaleprefs)))
  datam$matemv_randpref<-10*(-1*(datam$matemv_randpref)+maxeuc)/maxeuc
  
  datam$matemv_randtrait<-apply(datamrand[,25:29],1,function(x) dist(rbind(x,maleprefs)))
  datam$matemv_randtrait<-10*(-1*(datam$matemv_randtrait)+maxeuc)/maxeuc
  
  
  #Females#
  #Compute self mate value for females
  dataf$selfmv_main<-apply(dataf[,20:24],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_main<-10*(-1*(dataf$selfmv_main)+maxeuc)/maxeuc
  
  #Compute self mate value based on a random preference point
  dataf$selfmv_randpref<-apply(dataf[,20:24],1,function(x) dist(rbind(x,randmaleprefs)))
  dataf$selfmv_randpref<-10*(-1*(dataf$selfmv_randpref)+maxeuc)/maxeuc
  
  #Compute self mate value using scrambled traits
  dataf$selfmv_randtrait<-apply(datafrand[,20:24],1,function(x) dist(rbind(x,maleprefs)))
  dataf$selfmv_randtrait<-10*(-1*(dataf$selfmv_randtrait)+maxeuc)/maxeuc
  
  #Compute partner mate value for females
  dataf$matemv_main<-apply(dataf[,25:29],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_main<-10*(-1*(dataf$matemv_main)+maxeuc)/maxeuc
  
  #Compute partner mate value based on random preferences
  dataf$matemv_randpref<-apply(dataf[,25:29],1,function(x) dist(rbind(x,randfemaleprefs)))
  dataf$matemv_randpref<-10*(-1*(dataf$matemv_randpref)+maxeuc)/maxeuc
  
  dataf$matemv_randtrait<-apply(datafrand[,25:29],1,function(x) dist(rbind(x,femaleprefs)))
  dataf$matemv_randtrait<-10*(-1*(dataf$matemv_randtrait)+maxeuc)/maxeuc
  
  
  
  ###Desirabilities###
  #Create dataframes of desirabilities for participants
  selfdesm_main<-data.frame(t(apply(datam[,20:24],1,function(x) 6-abs(x-femaleprefs))))
  selfdesm_randpref<-data.frame(t(apply(datam[,20:24],1,function(x) 6-abs(x-randfemaleprefs))))
  selfdesm_randtrait<-data.frame(t(apply(datamrand[,20:24],1,function(x) 6-abs(x-femaleprefs))))
  
  matedesm_main<-data.frame(t(apply(datam[,25:29],1,function(x) 6-abs(x-maleprefs))))
  matedesm_randpref<-data.frame(t(apply(datam[,25:29],1,function(x) 6-abs(x-randmaleprefs))))
  matedesm_randtrait<-data.frame(t(apply(datamrand[,25:29],1,function(x) 6-abs(x-maleprefs))))
  
  selfdesf_main<-data.frame(t(apply(dataf[,20:24],1,function(x) 6-abs(x-maleprefs))))
  selfdesf_randpref<-data.frame(t(apply(dataf[,20:24],1,function(x) 6-abs(x-randmaleprefs))))
  selfdesf_randtrait<-data.frame(t(apply(datafrand[,20:24],1,function(x) 6-abs(x-maleprefs))))
  
  matedesf_main<-data.frame(t(apply(dataf[,25:29],1,function(x) 6-abs(x-femaleprefs))))
  matedesf_randpref<-data.frame(t(apply(dataf[,25:29],1,function(x) 6-abs(x-randfemaleprefs))))
  matedesf_randtrait<-data.frame(t(apply(datafrand[,25:29],1,function(x) 6-abs(x-femaleprefs))))
  
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
  
  #Proportion of variance explained by d-factors for males
  fvar_main<-ffa_main$Vaccounted[2]
  fvar_randpref<-ffa_randpref$Vaccounted[2]
  fvar_randtrait<-ffa_randtrait$Vaccounted[2]
  
  
  #Average variance explained across males and females
  var_main<-(mvar_main+fvar_main)/2
  var_randpref<-(mvar_randpref+fvar_randpref)/2
  var_randtrait<-(mvar_randtrait+fvar_randtrait)/2
  
  #Save to the appropriate row in the vardata dataframe
  vardata$fa_main[i]<-var_main
  vardata$fa_randpref[i]<-var_randpref
  vardata$fa_randtrait[i]<-var_randtrait
  
  
  ###d-Factor Loading as a Function of Cross-Character Assortment###
  
  #Calculate how much each trait generates cross-character assortment for males
  mcca_main<-sapply(1:ncol(selfdesm_main),function(x) mean(cor(selfdesm_main[,x],matedesm_main[,-x],use="pairwise.complete.obs")))
  mcca_randpref<-sapply(1:ncol(selfdesm_randpref),function(x) mean(cor(selfdesm_randpref[,x],matedesm_randpref[,-x],use="pairwise.complete.obs")))
  mcca_randtrait<-sapply(1:ncol(selfdesm_randtrait),function(x) mean(cor(selfdesm_randtrait[,x],matedesm_randtrait[,-x],use="pairwise.complete.obs")))
  
  #Calculate how much each trait generates cross-character assortment for females
  fcca_main<-sapply(1:ncol(selfdesf_main),function(x) mean(cor(selfdesf_main[,x],matedesf_main[,-x],use="pairwise.complete.obs")))
  fcca_randpref<-sapply(1:ncol(selfdesf_randpref),function(x) mean(cor(selfdesf_randpref[,x],matedesf_randpref[,-x],use="pairwise.complete.obs")))
  fcca_randtrait<-sapply(1:ncol(selfdesf_randtrait),function(x) mean(cor(selfdesf_randtrait[,x],matedesf_randtrait[,-x],use="pairwise.complete.obs")))
  
  #Save cross-character assortment values to the cca dataframe
  ccadata$cca[ccadata$CIN==i]<-c(fcca_main,mcca_main,fcca_randpref,mcca_randpref,fcca_randtrait,mcca_randtrait)
  
  #Extract the factor loading for each trait
  ccadata$loading[ccadata$CIN==i]<-c(as.numeric(ffa_main$loading),as.numeric(mfa_main$loading),as.numeric(ffa_randpref$loading),as.numeric(mfa_randpref$loading),as.numeric(ffa_randtrait$loading),as.numeric(mfa_randtrait$loading))
  
  
  
  ###Save Factor Scores###
  datam$fa_main<-as.numeric(mfa_main$scores)
  datam$fa_randpref<-as.numeric(mfa_randpref$scores)
  datam$fa_randtrait<-as.numeric(mfa_randtrait$scores)
  
  dataf$fa_main<-as.numeric(ffa_main$scores)
  dataf$fa_randpref<-as.numeric(ffa_randpref$scores)
  dataf$fa_randtrait<-as.numeric(ffa_randtrait$scores)
  
  datam<-cbind(datam,selfdesm_main)
  dataf<-cbind(dataf,selfdesf_main)
  
  colnames(datam[,39:43])<-c("des_intelligence","des_kindness","des_health","des_physatt","des_resources")
  colnames(dataf[,39:43])<-c("des_intelligence","des_kindness","des_health","des_physatt","des_resources")
  
  #Recombine the data
  loopdata<-rbind(dataf,datam)
  
  #Save the data
  procdata[procdata$CIN==i,]<-loopdata
  
}

#Rename procdata
data<-procdata

#Make country ID a factor
data$CIN<-as.factor(data$CIN)

#Make country ID a factor
ccadata$CIN<-as.factor(ccadata$CIN)


###Data Analysis###

##Trait and Desirability Correlations##

#Calculate correlations between traits and desirabilites separately for males and females
mcors<-cor(data[data$sex==1,c(20:24,39:43)],use="pairwise.complete.obs")
fcors<-cor(data[data$sex==0,c(20:24,39:43)],use="pairwise.complete.obs")

#Compile these into a single correlation matrix, with females in the upper triangle
cormatrix<-matrix(NA,10,10)
cormatrix[lower.tri(cormatrix)]<-mcors[lower.tri(mcors)]
cormatrix[upper.tri(cormatrix)]<-fcors[upper.tri(fcors)]

colnames(cormatrix)<-colnames(fcors)
rownames(cormatrix)<-rownames(fcors)

##Factor Variance Explained##

#Average variance explained across cultures
var<-data.frame("factor"=c("fa_main","fa_randpref","fa_randtrait"))
var$varexp<-c(mean(vardata$fa_main),mean(vardata$fa_randpref),mean(vardata$fa_randtrait))
var$lci<-var$varexp-c((sd(vardata$fa_main)/sqrt(nrow(vardata)))*qt(.975,nrow(vardata)-1),(sd(vardata$fa_randpref)/sqrt(nrow(vardata)))*qt(.975,nrow(vardata)-1),(sd(vardata$fa_randtrait)/sqrt(nrow(vardata)))*qt(.975,nrow(vardata)-1))
var$uci<-var$varexp+c((sd(vardata$fa_main)/sqrt(nrow(vardata)))*qt(.975,nrow(vardata)-1),(sd(vardata$fa_randpref)/sqrt(nrow(vardata)))*qt(.975,nrow(vardata)-1),(sd(vardata$fa_randtrait)/sqrt(nrow(vardata)))*qt(.975,nrow(vardata)-1))


##Factor Scores and Mate Value##
facmv_main<-lmer(fa_main~selfmv_main+(1+selfmv_main|CIN),data=data,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))
facmv_randpref<-lmer(fa_randpref~selfmv_randpref+(1+selfmv_randpref|CIN),data=data,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))
facmv_randtrait<-lmer(fa_randtrait~selfmv_randtrait+(1+selfmv_randtrait|CIN),data=data,control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

facmvfit_main<-predict(facmv_main)
facmvfit_randpref<-predict(facmv_randpref)
facmvfit_randtrait<-predict(facmv_randtrait)

mvslopes_main<-coef(facmv_main)$CIN
mvslopes_randpref<-coef(facmv_randpref)$CIN
mvslopes_randtrait<-coef(facmv_randtrait)$CIN

mvslopes_main<-mvslopes_main[order(mvslopes_main$selfmv_main),]
mvslopes_randpref<-mvslopes_randpref[order(mvslopes_randpref$selfmv_randpref),]
mvslopes_randtrait<-mvslopes_randtrait[order(mvslopes_randtrait$selfmv_randtrait),]

facmvcoef_main<-c(summary(facmv_main)$coefficients[1,1],summary(facmv_main)$coefficients[2,1])
facmvcoef_randpref<-c(summary(facmv_randpref)$coefficients[1,1],summary(facmv_randpref)$coefficients[2,1])
facmvcoef_randtrait<-c(summary(facmv_randtrait)$coefficients[1,1],summary(facmv_randtrait)$coefficients[2,1])

#Standardized slopes
#If Satterthwaite's approximation fails, use anova(model,ddf="Kenward-Roger")
facmv_mainstd<-stdcoef(facmv_main)
facmv_randprefstd<-stdcoef(facmv_randpref)
facmv_randtraitstd<-stdcoef(facmv_randtrait)



##Desirability Predictive Power and Factor Loading
cca_main<-lmer(loading~cca+(1+cca|CIN),data=ccadata[ccadata$factor=="fa_main",],control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))
cca_randpref<-lmer(loading~cca+(1+cca|CIN),data=ccadata[ccadata$factor=="fa_randpref",],control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))
cca_randtrait<-lmer(loading~cca+(1+cca|CIN),data=ccadata[ccadata$factor=="fa_randtrait",],control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))

#Find the overall line of best fit across model runs
ccaovlfit_main<-predict(cca_main)
ccaovlfit_randpref<-predict(cca_randpref)
ccaovlfit_randtrait<-predict(cca_randtrait)

#Find the coefficients for each model run's regression line
ccacoef_main<-c(summary(cca_main)$coefficients[1,1],summary(cca_main)$coefficients[2,1])
ccacoef_randpref<-c(summary(cca_randpref)$coefficients[1,1],summary(cca_randpref)$coefficients[2,1])
ccacoef_randtrait<-c(summary(cca_randtrait)$coefficients[1,1],summary(cca_randtrait)$coefficients[2,1])

#Standardized coefficients
cca_mainstd<-stdcoef(cca_main)
cca_randprefstd<-stdcoef(cca_randpref)
cca_randtraitstd<-stdcoef(cca_randtrait)



###Plots###

###Variance Explained by d-Factor###
#Create the plotting dataframe
vardata$country<-unique(data$country)
varplotdata<-melt(vardata,id.vars=c("CIN","country"))

varyaxis<-expression(paste("Variance Explaned By",italic(" d"),"-Factor"))

#Plots
varplot_main<-qplot(country,value,data=varplotdata[varplotdata$variable=="fa_main",],fill=country,xlab="Country",ylab=varyaxis,geom="blank")+geom_bar(stat="identity",position="dodge")+theme_classic(base_size=12)+theme(legend.position="none")+scale_fill_discrete(name="Factor")+geom_hline(yintercept=var[1,2],size=2)+theme(axis.text.x=element_text(angle=90,vjust=0,hjust=1))+annotate("text",x=2,y=max(varplotdata[varplotdata$variable=="fa_main",]$value)*1.05,fontface="bold",size=8,label="A")
varplot_randpref<-qplot(country,value,data=varplotdata[varplotdata$variable=="fa_randpref",],fill=country,xlab="Country",ylab=varyaxis,geom="blank")+geom_bar(stat="identity",position="dodge")+theme_classic(base_size=12)+theme(legend.position="none")+scale_fill_discrete(name="Factor")+geom_hline(yintercept=var[2,2],size=2)+theme(axis.text.x=element_text(angle=90,vjust=0,hjust=1))+annotate("text",x=2,y=max(varplotdata[varplotdata$variable=="fa_randpref",]$value)*1.05,fontface="bold",size=8,label="A")
varplot_randtrait<-qplot(country,value,data=varplotdata[varplotdata$variable=="fa_randtrait",],fill=country,xlab="Country",ylab=varyaxis,geom="blank")+geom_bar(stat="identity",position="dodge")+theme_classic(base_size=12)+theme(legend.position="none")+scale_fill_discrete(name="Factor")+geom_hline(yintercept=var[3,2],size=2)+theme(axis.text.x=element_text(angle=90,vjust=0,hjust=1))+annotate("text",x=2,y=max(varplotdata[varplotdata$variable=="fa_randtrait",]$value)*1.05,fontface="bold",size=8,label="A")


###d-Factor Score as a Function of Mate Value###
#Create plotting dataframes
mvplotdata_main<-data[complete.cases(data[,c(2,30,36)]),c(2,30,36)]
mvplotdata_randpref<-data[complete.cases(data[,c(2,31,37)]),c(2,31,37)]
mvplotdata_randtrait<-data[complete.cases(data[,c(2,32,38)]),c(2,32,38)]

mvplotyaxis<-expression(paste(italic("d"),"-Factor Score"))

#Plots
facmvplot_main<-qplot(selfmv_main,fa_main,color=as.factor(CIN),data=mvplotdata_main,xlab="Participant Mate Value",ylab=mvplotyaxis)+geom_line(aes(y=facmvfit_main))+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(mvplotdata_main$selfmv_main)*.95,y=max(mvplotdata_main$fa_main)*.9,fontface="bold",size=8,label="B")+geom_abline(intercept=facmvcoef_main[1],slope=facmvcoef_main[2],size=2.5,color="black")
facmvplot_randpref<-qplot(selfmv_randpref,fa_randpref,color=as.factor(CIN),data=mvplotdata_randpref,xlab="Participant Mate Value",ylab=mvplotyaxis)+geom_line(aes(y=facmvfit_randpref))+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(mvplotdata_randpref$selfmv_randpref)*.95,y=max(mvplotdata_randpref$fa_randpref)*.9,fontface="bold",size=8,label="B")+geom_abline(intercept=facmvcoef_randpref[1],slope=facmvcoef_randpref[2],size=2.5,color="black")
facmvplot_randtrait<-qplot(selfmv_randtrait,fa_randtrait,color=as.factor(CIN),data=mvplotdata_randtrait,xlab="Participant Mate Value",ylab=mvplotyaxis)+geom_line(aes(y=facmvfit_randtrait))+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(mvplotdata_randtrait$selfmv_randtrait)*.95,y=max(mvplotdata_randtrait$fa_randtrait)*.9,fontface="bold",size=8,label="B")+geom_abline(intercept=facmvcoef_randtrait[1],slope=facmvcoef_randtrait[2],size=2.5,color="black")


###d-Factor Loading as a Function of Cross-Character Assortment###
#Plots

predyaxis<-expression(paste(italic("d"),"-Factor Loadings"))

ccaplot_main<-qplot(cca,loading,color=as.factor(CIN),data=ccadata[ccadata$factor=="fa_main",],xlab="Degree of Cross-Character Assortment",ylab=predyaxis)+geom_line(aes(y=ccaovlfit_main),size=1)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(ccadata[ccadata$factor=="fa_main",]$cca)*.95,y=max(ccadata[ccadata$factor=="fa_main",]$loading)*.95,fontface="bold",size=8,label="C")+geom_abline(intercept=ccacoef_main[1],slope=ccacoef_main[2],size=2.5,color="black")
ccaplot_randpref<-qplot(cca,loading,color=as.factor(CIN),data=ccadata[ccadata$factor=="fa_randpref",],xlab="Degree of Cross-Character Assortment",ylab=predyaxis)+geom_line(aes(y=ccaovlfit_randpref),size=1)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(ccadata[ccadata$factor=="fa_randpref",]$cca)*1.1,y=max(ccadata[ccadata$factor=="fa_randpref",]$loading)*1.05,fontface="bold",size=8,label="C")+geom_abline(intercept=ccacoef_randpref[1],slope=ccacoef_randpref[2],size=2.5,color="black")
ccaplot_randtrait<-qplot(cca,loading,color=as.factor(CIN),data=ccadata[ccadata$factor=="fa_randtrait",],xlab="Degree of Cross-Character Assortment",ylab=predyaxis)+geom_line(aes(y=ccaovlfit_randtrait),size=1)+theme_classic(base_size=12)+theme(legend.position="none")+annotate("text",x=min(ccadata[ccadata$factor=="fa_randtrait",]$cca)*1.1,y=max(ccadata[ccadata$factor=="fa_randtrait",]$loading)*.95,fontface="bold",size=8,label="C")+geom_abline(intercept=ccacoef_randtrait[1],slope=ccacoef_randtrait[2],size=2.5,color="black")


plot_main<-grid.arrange(varplot_main,facmvplot_main,ccaplot_main,ncol=2,layout_matrix=cbind(c(1,1),c(2,3)),widths=c(1,.75))
plot_randpref<-grid.arrange(varplot_randpref,facmvplot_randpref,ccaplot_randpref,ncol=2,layout_matrix=cbind(c(1,1),c(2,3)),widths=c(1,.75))
plot_randtrait<-grid.arrange(varplot_randtrait,facmvplot_randtrait,ccaplot_randtrait,ncol=2,layout_matrix=cbind(c(1,1),c(2,3)),widths=c(1,.75))
cc_loadinghist<-qplot(ccadata$loading[ccadata$factor=="fa_main"],xlab="Desirability Factor Loading",ylab="Frequency",fill=I("coral"),color=I("black"),binwidth=.05)+theme_classic(base_size=20)

