rm(list=ls())
library("lme4")
library("multcomp")
library("lattice")
library(nlme)
library(survival)

U1611=read.csv("U1611_full.csv",h=T)
days = c(9,12,15,19,22,26,29,33,36)
U1611.long = reshape(idvar='Temp.ID',varying = c(7:15),timevar = 'Day',times=days,data=U1611,direction = 'long',v.names ='TV')
U1611.trt = read.csv("U1611_treatment.csv",h=T)
U1611.full = merge(U1611.long,U1611.trt)
U1611=merge(U1611,U1611.trt)
U1611.full = U1611.full[order(U1611.full$Group,U1611.full$Temp.ID,U1611.full$Day),]
#U1611.full$Treatment = as.factor(paste(U1611.full$Group,U1611.full$Treatment,sep=':'))
U1611.full$Treatment = as.factor(U1611.full$Treatment)
#U1611.full$Treatment = factor(U1611.full$Treatment,order=T,levels=levels(U1611.full$Treatment)[c(2,4,3,5,6,1)])
U1611.full$logTV = log(U1611.full$TV)
names(U1611.full)[2]="mouse"
attach(U1611.full)
tapply(TV,list(Treatment,Day),mean,na.rm=T)
boxplot(TV~Day) #approximately normal at specific day, 
plot(log(TV)~Day) #log(TV) linearly correlated with Day

U1611.gr<-groupedData(TV~Day | mouse, outer=~Group, order.groups=F,data=U1611.full)
plot(U1611.gr, display='mouse', outer=TRUE, aspect=2, key=F, xlab='Day', ylab='TV', main='Individual TV data by Treatment Group')
U1611.glog<-groupedData(logTV~Day | mouse, outer=~Group, order.groups=F,data=U1611.full)
plot(U1611.glog, display='mouse', outer=TRUE, aspect=2, key=F, xlab='Day', ylab='TV', main='Individual TV data by Treatment Group')


U1611.full$Group = as.factor(U1611.full$Group)
U1611.full$Anti.PD.L1 = as.factor(U1611.full$Anti.PD.L1)
U1611.full$Anti.PD.1 = as.factor(U1611.full$Anti.PD.1)
d26 = subset(U1611.full,Day==26)
#d26$Group = as.factor(d26$Group)
#d26$Anti.PD.L1 = as.factor(d26$Anti.PD.L1)
#d26$Anti.PD.1 = as.factor(d26$Anti.PD.1)
#m_aov = aov(TV~Group,data=d26)
#summary(m_aov)
#TukeyHSD(m_aov)
m_aov_log = aov(logTV~Group,data=d26)
summary(m_aov_log)
TukeyHSD(m_aov_log)

#RM anova
a2 = function(d){
	print(leveneTest(logTV~as.factor(Anti.PD.L1)*Treatment,data=d))
	aov.2=with(d,aov(logTV~Treatment*Anti.PD.L1))
	print(tapply(d[,"logTV"],list(d[,"Treatment"],d[,"Anti.PD.L1"]),mean,na.rm=T))
	print(summary(aov.2))
	print(TukeyHSD(aov.2))
	par(mar=c(5.1,20,3.1,1))
	plot(TukeyHSD(aov.2),las=1)
	par(mar=c(5.1,5.1,3.1,1))
}
g1 = subset(d26,Group %in% 1:6) #Group 1
a2(g1)
with(g1,pairwise.t.test(logTV,Group,p.adjust.method="bonferroni"))
model.tables(g1.aov.2,type="means",se=T)
require(gplots)
plotmeans(g1$TV~g1$Group,xlab="Group",ylab="TV", main="Mean Plot\nwith 95% CI")
g2 = subset(d26,Group %in% c(1,4,7,8))
a2(g2)
a3 = function(d){
	print(leveneTest(logTV~as.factor(Anti.PD.1)*Treatment,data=d))
	aov.2=with(d,aov(logTV~Treatment*Anti.PD.1))
	print(tapply(d[,"logTV"],list(d[,"Treatment"],d[,"Anti.PD.1"]),mean,na.rm=T))
	print(summary(aov.2))
	print(TukeyHSD(aov.2))
	par(mar=c(5.1,15,3.1,1))
	plot(TukeyHSD(aov.2),las=1)
	par(mar=c(5.1,5.1,3.1,1))
}
g3 = subset(d26,Group %in% c(1,9,14:17))
a3(g3)
g4 = subset(d26,Group %in% c(1,9,18:21))
a3(g4)
#RM ANOVA
rma = function(d){
	print(leveneTest(logTV~as.factor(Anti.PD.L1)*Treatment*as.factor(Day),data=d))
	aov.2=with(d,aov(logTV~Treatment*as.factor(Anti.PD.L1)*as.factor(Day)+Error(mouse/as.factor(Day))))
	print(summary(aov.2))
}
g1 = subset(U1611.full,Group %in% 1:6) #Group 1
rma(g1)
g2 = subset(U1611.full,Group %in% c(1,4,7,8)) #Group 1
rma(g2)
rma3 = function(d){
	print(leveneTest(logTV~as.factor(Anti.PD.1)*Treatment*as.factor(Day),data=d))
	aov.2=with(d,aov(logTV~Treatment*as.factor(Anti.PD.1)*as.factor(Day)+Error(mouse/as.factor(Day))))
	print(summary(aov.2))
	#print(model.tables(aov.2,type="means",se=T))
}
g3 = subset(U1611.full,Group %in% c(1,9,14:17))
rma3(g3)
g4 = subset(U1611.full,Group %in% c(1,9,18:21))
rma3(g4)


#Linear mixed model (to be tested)
U1611.noNA = U1611.full[complete.cases(U1611.full),]
g1 = subset(U1611.full,Group %in% 1:6) #Group 1

g1.gp<-groupedData(TV~Day | mouse, outer=~Treatment*Anti.PD.L1, order.groups=F,data=g1)
plot(g1.gp, display='mouse', outer=TRUE, aspect=2, key=F, xlab='Day', ylab='TV', main='Individual TV data by Treatment Group')

g1.gr<-groupedData(TV~Day|mouse, data=g1, order.groups=F)
m8<-lme(TV~Day*Treatment*Anti.PD.L1+I(Day^2)+I(Day^2):Treatment+I(Day^2):Anti.PD.L1, random=~Day-1, method='REML', data=g1.gr)
m9<-lme(TV~Day*Treatment*Anti.PD.L1+I(Day^2)+I(Day^2):Treatment+I(Day^2):Anti.PD.L1, random=~I(Day^2)+Day-1, method='REML', data=g1.gr)
#m10<-lme(TV~Day*Treatment*Anti.PD.L1+I(Day^2)+I(Day^2):Treatment+I(Day^2):Anti.PD.L1, random=~I(Day^2)+Day, method='REML', data=g1.gr)
summary(m8)
summary(m9)
anova(m8,m9)
anova(m9)
#when running lme, using different treatments as reference, the resulted p-value is different
plot(augPred(m9,level=0:1),layout=c(10,6,1),main="Individual profile & Mean profile") #model bulit using Raw TV

lmm = function(g){
	g.gr<-groupedData(TV~Day|mouse, data=g, order.groups=F)
	m8<-lme(TV~Day*Treatment*Anti.PD.L1+I(Day^2)+I(Day^2):Treatment+I(Day^2):Anti.PD.L1, random=~Day-1, method='REML', data=g.gr, na.action=na.omit)
	m9<-lme(TV~Day*Treatment*Anti.PD.L1+I(Day^2)+I(Day^2):Treatment+I(Day^2):Anti.PD.L1, random=~I(Day^2)+Day-1, method='REML', data=g.gr,na.action=na.omit)
	print(anova(m8,m9))
	print(anova(m9))
	print(summary(m9))
}
g2 = subset(U1611.full,Group %in% c(1,4,7,8)) #Group 2
lmm(g2)

lmm3 = function(g){
	g.gr<-groupedData(TV~Day|mouse, data=g, order.groups=F)
	m8<-lme(TV~Day*Treatment*Anti.PD.1+I(Day^2)+I(Day^2):Treatment+I(Day^2):Anti.PD.1, random=~Day-1, method='REML', data=g.gr, na.action=na.omit)
	m9<-lme(TV~Day*Treatment*Anti.PD.1+I(Day^2)+I(Day^2):Treatment+I(Day^2):Anti.PD.1, random=~I(Day^2)+Day-1, method='REML', data=g.gr,na.action=na.omit)
	print(anova(m8,m9))
	print(anova(m9))
	print(summary(m9))
}
g3 = subset(U1611.full,Group %in% c(1,9,14:17))#Group 2
lmm3(g3)
lmm4 = function(g){
	g.gr<-groupedData(TV~Day|mouse, data=g, order.groups=F)
	m8<-lme(TV~Day*Treatment*Anti.PD.1+I(Day^2)+I(Day^2):Treatment+I(Day^2):Anti.PD.1, random=~Day-1, method='REML', data=g.gr, na.action=na.omit)
	#m9<-lme(TV~Day*Treatment*Anti.PD.1+I(Day^2)+I(Day^2):Treatment+I(Day^2):Anti.PD.1, random=~Day, method='REML', data=g.gr,na.action=na.omit)
	#print(anova(m8,m9))
	print(anova(m8))
	print(summary(m8))
}
g4 = subset(U1611.full,Group %in% c(1,9,18:21))
lmm4(g4)

#survival analysis (not tested)
status = ifelse(U1611$Code=="TS" & U1611$X18 < 2000,0,1)
U1611$status = status
s1=with(Surv(Day,status),data=U1611)
survfit(s1~Treatment,data=U1611)
with(survdiff(s1~Treatment,rho=1),data=U1611)
with(survdiff(s1~Treatment,rho=0),data=U1611) #log rank test
s2=with(coxph(s1~Treatment),data=U1611)
summary(s2)
library(coin)
logrank_test(s1~Treatment,data=U1611)
U1611_12 = U1611[U1611$Group %in% c(1,2),]
U1611_13 = U1611[U1611$Group %in% c(1,3),]
U1611_14 = U1611[U1611$Group %in% c(1,4),]
U1611_15 = U1611[U1611$Group %in% c(1,5),]
U1611_16 = U1611[U1611$Group %in% c(1,6),]
logrank_test(Surv(Day,status)~Treatment,data=U1611_12,distribution = "exact")
logrank_test(Surv(Day,status)~Treatment,data=U1611_13,distribution = "exact")
logrank_test(Surv(Day,status)~Treatment,data=U1611_14,distribution = "exact")
logrank_test(Surv(Day,status)~Treatment,data=U1611_15,distribution = "exact")
logrank_test(Surv(Day,status)~Treatment,data=U1611_16,distribution = "exact")



