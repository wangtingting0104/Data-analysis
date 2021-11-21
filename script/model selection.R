rm(list=ls())
library(metafor)
library(glmulti)

gcdata<-read.csv("RESPONSE.csv",header=T,row.names = 1,check.names=F)

para<-c("richness")##richness,shannon,evenness,PD,Nitrification,nitrogen_fixation,sulfur,cellulolysis,denitrification,chitinolysis,xylanolysis,ureolysis,methol,ligninolysis,metholotrophy,pathogen,hetero,toxin,sporulation
i=1

mydata<-gcdata[,c("pHRR","TNRR","SOCRR","Ecosystems","TFamily","Longitude","Latitude",
                  para[i],paste(para[i],"_var",sep=""))]#
mydata<-na.omit(mydata)
names(mydata)[ncol(mydata)-1]<-c("yi")
names(mydata)[ncol(mydata)]<-c("vi")
##rma
##
rma.glmulti <- function(formula, data, ...) { rma(formula, vi, data=data, method="ML", ...)} ##
##
res <- glmulti(yi ~ pHRR+TNRR+SOCRR+Ecosystems+TFamily+Longitude+Latitude, data=mydata,
               level=1, fitfunction=rma.glmulti, crit="aicc", confsetsize=128)#confsetsize

plot(res)  # AICC of different model##
plot(res, type="s")  # Importance of different predictors##

reswei <- weightable(res)##
reswei
summary(res@objects[[1]])# show the detail of specific model of order of 1##
transmodel <- reswei[reswei$aicc <= min(reswei$aicc) + 2,]  # best model (difference in aicc within 2)##输出最优模型，有时候最优模型不止一个，可以有多个模型与AICc值与AICc值最低的模型的AICc差值<2

eval(metafor:::.glmulti)
transres <- as.data.frame(coef(res))  # Weighted averages of the model coefficients
transres$nobser <- nrow(mydata)
transres$nfactor <- ncol(mydata)-2
transres$vari <- para[i]
transres$factor<-row.names(transres)
#transmodel$nobser <- nrow(mydata)
#transmodel$nfactor <- ncol(mydata)-2
#transmodel$vari <- para[i]
modelres<-transres#
#bestmodel<-transmodel#

write.csv(modelres, paste(para,"modelres.csv"), row.names = F, na="")  
#write.csv(bestmodel, paste(para,"bestmodel.csv"), row.names = F, na="")

