run=function(ichunk){
	print(ichunk)
	#rm(list=ls())
	library(frailtypack)
	library(frailtySurv)
	library(snow)
	simuFrailty = function(n_pdx=10,n_mouse=5,b=-0.693,theta=1.538,scale=0.0545,shape=3.72,censor_day=21,ns=100){
		m=rep(rep(c(0,1),c(n_mouse,n_mouse)),n_pdx)
		nZ = rep(0,ns)
		cut_off = qnorm(0.025)
		for(rep in 1:ns){
			sim_data=genfrail(N = n_pdx, K = 2*n_mouse,beta=b,frailty='gamma',theta=theta, Lambda_0=function(t, tau=shape, C=scale) (C*t)^tau,covar.matrix=as.matrix(m))
			sim_data$Day = sim_data$time
			sim_data$status = ifelse(sim_data$Day<=censor_day,1,0)
			names(sim_data) = c('PDX','Mouse','time','censoring','Treatment','Day')
			tryCatch(f1<-frailtyPenal(Surv(Day,censoring)~cluster(PDX)+Treatment,data=sim_data,hazard="Weibull"),error=function(e) NULL )
			f1 = get0("f1",ifnotfound=NULL)
			if(is.null(f1)){
				#print("Oops!",n_pdx,n_mouse,b)
				nZ[rep]=Inf
			}else{
				nZ[rep]=f1$coef/sqrt(f1$varH)
			}	
		}
		nZ=nZ[nZ!=-Inf & nZ != Inf]
		sum(nZ<cut_off)/ns
	}
	sample = list(n_pdx=2:3,n_mouse=1:2)
	grid = expand.grid(sample)
	n_grid = dim(grid)[1]
	out = list()
	for(bi in ichunk){
		results = matrix(rep(-1,4*n_grid),n_grid,4)
		colnames(results)=c('n_pdx','n_mouse','bt','Power')
		fn = paste('simu_',as.character(bi),'.csv',sep='')
		write(colnames(results),fn,ncolumns=4,sep=',')
		for(i in 1:n_grid){
			s = grid[i,]
			powerE = simuFrailty(n_pdx=s$n_pdx,n_mouse=s$n_mouse,b=bi,theta=1,censor_day=21,ns=100)
			results[i,] = c(s$n_pdx,s$n_mouse,bi,powerE)
			write(results[i,],fn,append=T,ncolumns=4,sep=',')
		}
		out[[as.character(bi)]] = results
	}
	out
}

paraRun = function(cls,b_vec){
	n=length(b_vec)
	nc=length(cls)
	options(warn=-1)
	ichunks=split(b_vec,1:nc)
	print(ichunks)
	options(warn=0)
	result=clusterApply(cls,ichunks,run)
	print(result)
}
nc=3
b_vec = rev(log(seq(0.7,0.9,0.1)))
cl=makeCluster(type="SOCK",rep("localhost",nc))
paraRun(cl,b_vec)




