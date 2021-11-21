library(meta)

library(metafor)


d1<-read.csv("Response.csv")

d2<-escalc(measure="ROM",data=d1,m1i=RMEAN,sd1i=RSD,n1i=RN,m2i=BMEAN,sd2i=BSD,n2i=BN)


r11<-rma(yi,vi, data=d2, method="DL")
funnel(r11,main="For birds")
#summary(r11)


#regtest(r11, model="rma")
regtest(r11, model="rma") ##Egger's regression test for funnel plot asymmetry,for the meta-analytic models
##Egger's regression test for funnel plot asymmetry,for the meta-analytic models
#regtest(r22, model="rma")


###Trim and fill method
###Trim and fill: A simple funnel-plot-based method 
###of testing and adjusting for publication bias 
###in meta-analysis. Biometrics, 56(2), 455–463. 

summary(r11)
t<-trimfill(r11,estimator="R0")
funnel(t)

summary(t)

