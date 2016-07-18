/* 
INPUT PARAMETERS
nreps = number of simulations
nclst = number of clusters
numobs = number of observations per cluster
intvar1 = variance of random cluster effects (intercepts)
sigeps1 = residual variance 
*/
#n_pdx: The number of PDX model
#n_mouse: The number of mouses in each PDX model, assume balance design
#n_rep: The number of time points for each mouse, assuming the same for each mouse, starting from 0 with a step of 3
#Three level longtitudal experiments. Level 1: PDX,Level 2: Mouse(Covariate:Treatment),Level 3:Time(Day) Covariate(Tumor Volume)
#NOT USED: Model: TV/log(TV) = b0+b1*Treatment+b2*Day+b3*Treatment*Day+u0(j)+u1(j)*Day+u0(i|j)+u1(i|j)*Day+e(tij)
#Model: TV/log(TV) = b0+b1*Day+bt*Treatment*Day+u0(j)+u1(j)*Day+u0(i|j)+u1(i|j)*Day+e(tij)
#Assume covariance matrix D1(PDX), D2(Mouse) (independent) ; simple variance structure for error term (Equal and independent)
#We cared about b3=Treatment interaction with time
#u0,u1,rho1 : u0,u1->sd of random effect of PDX on intercept and Day, rho1->correlation
#r0,r1,rho2 : r0,r1->sd of random effect of mouse on intercept and Day, rho2->correlation
#eps: sd of error term
rm(list=ls())
library(MASS)
library(nlme)
simuLMM = function(n_pdx=5,n_mouse=5,n_rep=7,bt=-0.026521,u0=0.1416,u1=0.016558,rho1=0.163,r0=0.31908,r1=0.019781,rho2=0.079,eps=0.19645,b0=5.26441,b1=0.0706,ns=100,alpha=0.05){
	#print(paste("n_pdx=",n_pdx,",n_mouse=",n_mouse,",n_day=",n_rep,",bt=",bt))
	m_pdx = matrix(c(u0^2,rho1*u0*u1,rho1*u0*u1,u1^2),2,2) #Random effect matrix of PDX
	m_mouse = matrix(c(r0^2,rho2*r0*r1,rho2*r0*r1,r1^2),2,2) #Random effect matrix of mouse
	s = list(day=seq(0,by=3,length.out=n_rep),mouse=1:(n_mouse*2),pdx=1:n_pdx) #Assume we have two groups, one treatment group and one control
	data = expand.grid(s)
	data$mouse = data$mouse + n_mouse*(data$pdx-1) #assign each mouse a different number
	treatment_pdx = rep(c(0,1),c(n_mouse*n_rep,n_mouse*n_rep)) #0: vehicle, 1:treatment
	treatment = rep(treatment_pdx,n_pdx)
	data$treatment = treatment
	cntrl = lmeControl(maxIter=0,msMaxIter=0,niterEM=0,returnObject=T,opt="optim")
	n_obs = n_pdx*n_mouse*2*n_rep
	tv_s = matrix(rep(0,n_obs*ns),n_obs,ns) #ns is number of simulations
	for(rep in 1:ns){
		u_pdx = mvrnorm(n_pdx,c(0,0),m_pdx) #random sample of pdx random effect
		r_pdx = mvrnorm(2*n_pdx*n_mouse,c(0,0),m_mouse) #random sample of mouse random effect
		error = rnorm(n_obs,mean=0,sd=eps)
		tv = rep(0,n_obs)
		for (i in 1:n_obs){
			obs_i = data[i,]
			#print(obs_i)
			#tv[i] = b0 + b1*obs_i$treatment+b2*obs_i$day+bt*obs_i$treatment*obs_i$day #fixed effect
			tv[i] = b0 +b1*obs_i$day+bt*obs_i$treatment*obs_i$day #fixed effect
			tv[i] = tv[i]+u_pdx[obs_i$pdx,]%*%c(1,obs_i$day)+r_pdx[obs_i$mouse,]%*%c(1,obs_i$day)+error[i] #random effect
		}
		tv_s[,rep] = tv
	}	
	#print(tv_s)
	data$tv = tv_s[,1]
	tryCatch(model.fit<-lme(tv ~ day + day:treatment, random = list(pdx = ~day, mouse = ~day), data=data, method='REML',control=cntrl),error=function(e) NULL )
	sims1 = anova(model.fit)
	simS = apply(tv_s[,-1],2,function(y){
							data$tv=y
							auxFit = tryCatch(update(model.fit,data=data),error=function(e) NULL )
							anova(auxFit)
						})
	simS[[ns]]=sims1
	#simS
	pstatE = sapply(simS,function(x) x["day:treatment","p-value"])
	powerE = sum(pstatE<alpha)/ns
	print(paste("n_pdx=",n_pdx,",n_mouse=",n_mouse,",n_day=",n_rep,",bt=",bt,",Power=",powerE))
	powerE
}

sample = list(n_pdx=2:50,n_mouse=1:10,n_rep=7:7,bt=seq(-0.0706/10,-0.0706/10*9,by=-0.0706/10)) #starting: treatment reduces rate by 10% to 90%
grid = expand.grid(sample)
n_grid = dim(grid)[1]
results = matrix(rep(-1,5*n_grid),n_grid,5)
colnames(results)=c('n_pdx','n_mouse','n_rep','bt','Power')
write(colnames(results),"result.txt",ncolumns=5,sep=',')
for(i in 1:n_grid){
	s = grid[i,]
	#print(system.time(simuLMM(s$n_pdx,s$n_mouse,s$n_rep,s$bt,ns=100)))
	powerE = simuLMM(s$n_pdx,s$n_mouse,s$n_rep,s$bt,ns=100,alpha=0.01)
	results[i,] = c(s$n_pdx,s$n_mouse,s$n_rep,s$bt,powerE)
	write(results[i,],"result.txt",append=T,ncolumns=5,sep=',')
}
