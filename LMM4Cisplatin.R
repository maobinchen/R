###########################################################################################
# Linear mixed models for clustered longitudinal data of cisplatin
# Note:
# 1. Baseline TV (TV0) is not used as a covariate, but another measure of the outcome
###########################################################################################
rm(list=ls())
library(nlme)
library(lattice)

###########################################################################################
#--read in data, draw growth curves
###########################################################################################
CisplatinDS<-read.table('CisplatinDataSet.txt', head=T)
logTV<-log(CisplatinDS$TV+1)
boxplot(CisplatinDS$ERCC1FPKM)
CisplatinDS<-data.frame(CisplatinDS, logTV)
CisplatinDS <- within(CisplatinDS, Treatment <- relevel(Treatment, ref = 'Vehicle'))
CisplatinDS.Vehicle<-subset(CisplatinDS,CisplatinDS$Treatment == 'Vehicle')
CisplatinDS.Cisplatin<-subset(CisplatinDS,CisplatinDS$Treatment == 'Cisplatin')
CisplatinDS.GA<-subset(CisplatinDS,CisplatinDS$CancerType == 'GA')
CisplatinDS.ES<-subset(CisplatinDS,CisplatinDS$CancerType == 'ES')
CisplatinDS.LU<-subset(CisplatinDS,CisplatinDS$CancerType == 'LU')
pdf(file='GrowthCurves.pdf')
xyplot(logTV~Day|PDX, groups=Mouse, type='l', data=CisplatinDS.Vehicle, main='Vehicle', as.table=TRUE, index.cond=list(unique(CisplatinDS.Vehicle$PDX)), ylim=c(4,8.5), ylab=expression(paste('log(TV (mm'^'3', '))', sep='')))
xyplot(logTV~Day|PDX, groups=Mouse, type='l', data=CisplatinDS.Cisplatin, main='Cisplatin', as.table=TRUE, index.cond=list(unique(CisplatinDS.Cisplatin$PDX)), ylim=c(4,8.5),ylab=expression(paste('log(TV (mm'^'3', '))', sep='')))
dev.off()

###########################################################################################
#--build models
#--step1: fit a model with a loaded 'mean' structure
###########################################################################################
model1.fit<-lme(logTV ~ Day + ERCC1FPKM + CancerType + Treatment + Day:ERCC1FPKM + Day:CancerType + Day:Treatment + Day:Treatment:ERCC1FPKM + Day:CancerType:Treatment, 
				random = list(PDX = ~Day, Mouse = ~1), data=CisplatinDS, method='REML')
summary(model1.fit)
intervals(model1.fit)
random.effects(model1.fit)
model2.fit<-lme(logTV ~ Day + Day:Treatment, 
				random = list(PDX = ~Day, Mouse = ~Day), data=CisplatinDS, method='ML')
model3.fit<-lme(logTV ~ Day + Treatment + Day:Treatment, 
				random = list(PDX = ~Day, Mouse = ~Day), data=CisplatinDS, method='ML')
anova(model2.fit,model3.fit)				
model2.fit<-lme(logTV ~ Day + Day:Treatment, 
				random = list(PDX = ~Day, Mouse = ~Day), data=CisplatinDS, method='REML')
				
				

###########################################################################################
#--step2: select a structure for the random effects (model1 vs model1A)
###########################################################################################
#--model1A is specified by omitting the random mouse effects from Model 1.
model1A.fit<-lme(logTV ~ Day + ERCC1FPKM + CancerType + Treatment + Day:ERCC1FPKM + Day:CancerType + Day:Treatment, 
				random = list(PDX = ~Day), data=CisplatinDS, method='REML')
summary(model1A.fit)
#--compare model1 and model1A: the p-value must be divided by 2, smaller AIC is better, here model1 is better
anova(model1.fit, model1A.fit)

###########################################################################################
#--step3: select a covariance structure for the residuals
###########################################################################################
model2A.fit<-lme(logTV ~ Day + ERCC1FPKM + CancerType + Treatment + Day:ERCC1FPKM + Day:CancerType + Day:Treatment, 
				random = list(PDX = ~Day, Mouse = ~1), 
				corr=corCompSymm(0.5, form=~1|PDX/Mouse),
				weights=varIdent(form=~1|Day),
				data=CisplatinDS, method='REML')
summary(model2A.fit) #does not converge

model2B.fit<-lme(logTV ~ Day + ERCC1FPKM + CancerType + Treatment + Day:ERCC1FPKM + Day:CancerType + Day:Treatment, 
				random = list(PDX = ~Day, Mouse = ~1), 
				corr=corCompSymm(0.5, form=~1|PDX/Mouse),
				data=CisplatinDS, method='REML')
summary(model2B.fit) #converge
intervals(model2B.fit)
anova(model1.fit, model2B.fit) #no differece, keep the loaded model1


model2C.fit<-lme(logTV ~ Day + ERCC1FPKM + CancerType + Treatment + Day:ERCC1FPKM + Day:CancerType + Day:Treatment, 
				random = list(PDX = ~Day, Mouse = ~1), 
				weights=varIdent(form=~1|Day),
				data=CisplatinDS, method='REML')
summary(model2C.fit) #does not converge

###########################################################################################
#--step4: reduce the model by removing nonsignificant fixed effects
###########################################################################################
#--first refit by ML
model1.ml.fit<-lme(logTV ~ Day + ERCC1FPKM + CancerType + Treatment + Day:ERCC1FPKM + Day:CancerType + Day:Treatment, 
				random = list(PDX = ~Day, Mouse = ~1), data=CisplatinDS, method='ML')
summary(model1.ml.fit)
intervals(model1.fit)

#--refit by removing non-significant fixed effect of the two-way interactions, also by ML
model3.ml.fit<-lme(logTV ~ Day + ERCC1FPKM + CancerType + Treatment + Day:CancerType + Day:Treatment, 
				random = list(PDX = ~Day, Mouse = ~1), data=CisplatinDS, method='ML')
summary(model3.ml.fit)
anova(model1.ml.fit,model3.ml.fit) #no difference, so we keep the simpler model model3

#--final fitting
model3.fit<-lme(logTV ~ Day + ERCC1FPKM + CancerType + Treatment + Day:CancerType + Day:Treatment, 
				random = list(PDX = ~Day, Mouse = ~1), data=CisplatinDS, method='REML')
summary(model3.fit)
intervals(model3.fit)
random.effects(model3.fit)


###########################################################################################
#model diagnostic
###########################################################################################
model.fit<-model3.fit

pdf(file='Diagnostics.pdf')
#--get standardized residuals
plot(model.fit, resid(., type="p")~fitted(.), xlab="Predicted value", ylab="Residual", abline=0, cex=0.2)

#--get qq plot
qqnorm(model.fit, ~resid(.)) #layout=c(2,1), aspect = 2, id = 0.05) 
qqnorm(model.fit, ~ranef(.,level=1), layout=c(1,2), aspect = 1, id = 0.05)
plot(model.fit, logTV~fitted(.)|factor(PDX), xlab="Predicted value", ylab="logTV", abline=c(0,1), cex=0.2)

dev.off()

###########################################################################################
#--refit after removing ES0214, ES0204 and LU1204
###########################################################################################
CisplatinDS2<-CisplatinDS[(CisplatinDS$PDX != 'ES0214' & CisplatinDS$PDX != 'LU1204' & CisplatinDS$PDX != 'ES0204'),]
model3A.fit<-lme(logTV ~ Day + ERCC1FPKM + CancerType + Treatment + Day:CancerType + Day:Treatment, 
				random = list(PDX = ~Day, Mouse = ~1), data=CisplatinDS2, method='REML')
summary(model3A.fit)
intervals(model3A.fit)
#random.effects(model3A.fit)

model.fit<-model3A.fit
pdf(file='Diagnostics2.pdf')
#--get standardized residuals
plot(model.fit, resid(., type="p")~fitted(.), xlab="Predicted value", ylab="Residual", abline=0, cex=0.2)

#--get qq plot
qqnorm(model.fit, ~resid(.)) #layout=c(2,1), aspect = 2, id = 0.05) 
qqnorm(model.fit, ~ranef(.,level=1), layout=c(1,2), aspect = 1, id = 0.05)
plot(model.fit, logTV~fitted(.)|factor(PDX), xlab="Predicted value", ylab="logTV", abline=c(0,1), cex=0.2)

dev.off()