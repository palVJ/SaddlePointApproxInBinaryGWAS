#Estimate type I errors conditional on the number of cases with a binary and continuous covariate:

args = commandArgs(TRUE)
arID = as.integer(args[1])
sims = as.numeric(args[2])
maf = as.double(args[3])
Cases = as.numeric(args[4])
Controls= as.integer(args[5])
N = Cases + Controls
#there are 4 covariates in total: intercept, X1, X2 and G
Ncovs = 4
ccratio = Cases/N
log.p = FALSE
Cutoff = 4

#Choose beta_0 such that the prevalence is around 0.01
beta_0 = -5.6

#integrand = function(x,c){(1/(1+exp(-c-x)))*(sqrt(1/(2*pi)))*exp(-0.5*x^2)}
#Prevalence = 0.5*integrate(integrand, lower = -Inf,upper = Inf, c = beta_0)$value + 0.5*integrate(integrand, lower = -Inf,upper = Inf, c = beta_0+1)$value
Prevalence = 0.0109492


#Phenotype values
pheno = c(rep(0,Controls),rep(1,Cases))

#Conditional probabilities of X1 given Y:

#Probability of X1 = 1 given Y1 = 1
#X1prob1_case = 0.5*integrate(integrand,lower = -Inf, upper = Inf, c = beta_0+1)$value/(Prevalence)
X1prob1_case = 0.7260159

#Probability of X1 = 1 given Y1 = 0
#X1prob1_control = (0.5-0.5*integrate(integrand,lower = -Inf, upper = Inf, c = beta_0+1)$value)/(1-Prevalence)
X1prob1_control = 0.4974979

#Log density of X2 given Y1 = 1:
#Log_Dens_X2_given_Y1 = function(x,beta_0,Prevalence){return(log((0.5*integrand(x,c = beta_0)+0.5*integrand(x,c = beta_0+1))/(Prevalence)))}
Log_Dens_X2_given_Y1 = function(x,beta_0,Prevalence){return(log((0.5*(1/(1+exp(-beta_0-x)))*(sqrt(1/(2*pi)))*exp(-0.5*x^2)+0.5*(1/(1+exp(-beta_0-1-x)))*(sqrt(1/(2*pi)))*exp(-0.5*x^2))/(Prevalence)))}

#The derivative of log density:
Log_Dens_X2_given_Y1_prime = function(x,beta_0,Prevalence){return(-x+(1/(1/(1+exp(-beta_0-x))+1/(1+exp(-beta_0-x-1))))*(exp(-beta_0-x)/(1+exp(-beta_0-x))^2+exp(-beta_0-1-x)/(1+exp(-beta_0-1-x))^2))}

#Log density of X2 given Y1 = 0:
#Log_Dens_X2_given_Y0 = function(x,beta_0,Prevalence){return(log((dnorm(x)- 0.5*integrand(x,c = beta_0)-0.5*integrand(x,c = beta_0+1))/(1-Prevalence)))}
Log_Dens_X2_given_Y0 = function(x,beta_0,Prevalence){return(log((dnorm(x)- 0.5*(1/(1+exp(-beta_0-x)))*(sqrt(1/(2*pi)))*exp(-0.5*x^2)-0.5*(1/(1+exp(-beta_0-1-x)))*(sqrt(1/(2*pi)))*exp(-0.5*x^2))/(1-Prevalence)))}

#The derivative of log density
Log_Dens_X2_given_Y0_prime = function(x,beta_0,Prevalence){return(-x+(1/(1-0.5/(1+exp(-beta_0-x))-0.5/(1+exp(-beta_0-x-1))))*(-0.5*exp(-beta_0-x)/(1+exp(-beta_0-x))^2-0.5*exp(-beta_0-1-x)/(1+exp(-beta_0-x-1))^2))}

samplerX2_given_Y0 = arscpp::ars(f = Log_Dens_X2_given_Y0,f_prime = Log_Dens_X2_given_Y0_prime,xlb = -Inf,xrb = Inf, x = c(-2,0,2),beta_0 = beta_0, Prevalence = Prevalence)
samplerX2_given_Y1 = arscpp::ars(f = Log_Dens_X2_given_Y1,f_prime = Log_Dens_X2_given_Y1_prime,xlb = -Inf,xrb = Inf, x = c(-2,0,2),beta_0 = beta_0, Prevalence = Prevalence)
source("~/DoubleSaddleInGWAS.R")
library(SPAtest)

#start for-loop. Preallocate matrix
pvals = matrix(data= NA, ncol = 10,nrow = 10000)
AddToRow = 1
colnames(pvals) = c("SPA_LR","SPA_BN","SPASCC_LR","SPASCC_BN","DSPASCC_LR","DSPASCC_BN","DSPASCCFAST1_LR","DSPASCCFAST1_BN","DSPASCCFAST2_LR","DSPASCCFAST2_BN")



SimsDone = sims

for(j in 1:sims){
 
  #simulate X1:
  X1 = c(rbinom(n=Controls,size = 1, prob = X1prob1_control),rbinom(n=Cases,size = 1, prob = X1prob1_case)) 

  #simulate X2:
  X2 = c(samplerX2_given_Y0$sample(Controls),samplerX2_given_Y1$sample(Cases))

  #create covariate matrix
  cov = as.matrix(data.frame(intercept = rep(1,N),X1 = X1, X2 = X2))
  #simulate genotype values
  G = rbinom(N,2,maf)

  if(sum(G)>0){	
  	#Find the MLE estimates under the null distribution (gamma = 0):
  	obj.null = ScoreTest_NULL_Model_fast(pheno = pheno, cov = cov)

  	muhat = obj.null$muhat
  	G_tilde = G - (obj.null$XXWX_inv)%*%((obj.null$XW)%*%G)
  	X = cbind(G_tilde,cov)
  	u = crossprod(as.matrix(X),as.matrix(pheno-muhat))
  	var1<-sum(muhat*(1-muhat)*G_tilde^2)
	#Observed score test statistic:
	u1 = u[1]
	#Standardised score test asymptotically normal distributed.
 	score = u1/sqrt(var1)

  	if(abs(score) > Cutoff){
                
                u1_min = -sum(G*muhat)
                u1_max = sum(G*(1-muhat))  

  		#Compute p-values using doubleSaddle assuming continuous density and double saddle with discrete distributions, first and second continuity corrections
		#out.uni1DSPA = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat,u1,u1_max,X,method = "DSPA")
		#out.uni1DSPAFCC = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat,u1,u1_max,X,method = "DSPAFCC")
		out.uni1DSPASCC = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat,u1,u1_max,X,method = "DSPASCC")
		m1 = sum(muhat*G_tilde)
		q = sum(G_tilde*pheno)

		#Find closest gridpoint to -u1 such that P(U1 <= -u1) = P(U1 <= uinv) for u1 positive and P(U1 >= -u1) = P(U1 >= uinv) for u1 negative:
                u1_inv = u1 - sign(u1)*ceiling(abs(2*u1))
		qinv = u1_inv + m1

		#Compute variables for approximate computations of CGFs
		Var_U_beta = obj.null$XWX
		G12 = which(G > 0)
		muhat_m = muhat[G12]
		cov_m = cov[G12,]
		Var_U_star = matrix(0,ncol = Ncovs,nrow = Ncovs)
		Var_U_star[2:Ncovs,2:Ncovs] = Var_U_beta - t(muhat_m*(1-muhat_m)*cov_m)%*%cov_m
		
		out.uni1DSPASCC_FAST = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat_m,u1,u1_max,X[G12,],method = "DSPASCC",Var_U_star = Var_U_star)
	
		out.uni1SPA<-SPAtest:::getroot_K1(0, mu=muhat, g=G_tilde, q=q)
		out.uni2SPA<-SPAtest:::getroot_K1(0, mu=muhat, g=G_tilde, q=qinv)
		
		#out.uni1SPAFCC = getrootSaddleContCorr_K1(0,muhat,q,u1,u1_max,G_tilde,method = "SPAFCC")
		out.uni1SPASCC = getrootSaddleContCorr_K1(0,muhat,q,u1,u1_max,G_tilde,method = "SPASCC")

		if(-u1 > u1_max || -u1 < u1_min){
        	#if minus observed score test statistic is outside range of density, no saddlepoint exists and we do one-sided p-value calculations
			#out.uni2DSPA = list(root = NA,convergence = 0)
			#out.uni2DSPAFCC = list(root = NA,convergence = 0)
			out.uni2DSPASCC = list(root = NA,convergence = 0)
			#out.uni2SPAFCC = list(root = NA,convergence = 0)
			out.uni2SPASCC = list(root = NA,convergence = 0)
			out.uni2DSPASCC_FAST = list(root = NA,convergence = 0) 
        	}else{
						
			#out.uni2DSPA = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat,u1_inv,u1_max,X,method = "DSPA")
 	        	#out.uni2DSPAFCC = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat,u1_inv,u1_max,X,method = "DSPAFCC")
                	out.uni2DSPASCC = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat,u1_inv,u1_max,X,method = "DSPASCC")
			#out.uni2SPAFCC = getrootSaddleContCorr_K1(0,muhat,qinv,u1_inv,u1_max,G_tilde,method = "SPAFCC")
                	out.uni2SPASCC = getrootSaddleContCorr_K1(0,muhat,qinv,u1_inv,u1_max,G_tilde,method = "SPASCC")
			out.uni2DSPASCC_FAST = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat_m,u1_inv,u1_max,X[G12,],method = "DSPASCC",Var_U_star = Var_U_star)
		}


		#p1DSPA = Get_DoubleSaddle_Prob(out.uni1DSPA$root,X,muhat,u1,log.p=FALSE,method = "DSPA",pval.comp = "both")		
		#p2DSPA = Get_DoubleSaddle_Prob(out.uni2DSPA$root,X,muhat,u1_inv,log.p=FALSE,method = "DSPA",pval.comp = "both")

		#p1DSPAFCC<-Get_DoubleSaddle_Prob(out.uni1DSPAFCC$root,X,muhat,u,log.p=FALSE,method = "DSPAFCC",pval.comp = "both")
        	#p2DSPAFCC <- Get_DoubleSaddle_Prob(out.uni2DSPAFCC$root,X,muhat,u1_inv,log.p=FALSE,method = "DSPAFCC",pval.comp = "both")

		p1SPA = Get_Saddle_Prob_SPA(that = out.uni1SPA$root,u1 = u1, muhat = muhat, G_tilde = G_tilde, q = q,log.p=FALSE,pval.comp = "both")
		p2SPA = Get_Saddle_Prob_SPA(that = out.uni2SPA$root,u1 = u1_inv, muhat = muhat, G_tilde = G_tilde, q = qinv,log.p=FALSE,pval.comp = "both")

		#p1SPAFCC = Get_Saddle_ProbContCor(out.uni1SPAFCC$root,G_tilde,muhat,u1,q,log.p=FALSE,method = "SPAFCC",pval.comp = "both")
		#p2SPAFCC = Get_Saddle_ProbContCor(out.uni2SPAFCC$root,G_tilde,muhat,u1_inv,qinv,log.p=FALSE,method = "SPAFCC",pval.comp = "both")

		p1SPASCC = Get_Saddle_ProbContCor(out.uni1SPASCC$root,G_tilde,muhat,u1,q,log.p=FALSE,method = "SPASCC",pval.comp = "both")
                p2SPASCC = Get_Saddle_ProbContCor(out.uni2SPASCC$root,G_tilde,muhat,u1_inv,qinv,log.p=FALSE,method = "SPASCC",pval.comp = "both")

        	p1DSPASCC<-Get_DoubleSaddle_Prob(out.uni1DSPASCC$root,X,muhat,u1,log.p=FALSE,method = "DSPASCC",pval.comp = "both")
        	p2DSPASCC<-Get_DoubleSaddle_Prob(out.uni2DSPASCC$root,X,muhat,u1_inv,log.p=FALSE,method = "DSPASCC",pval.comp = "both")
                    
		p1DSPASCC_FAST1<-Get_DoubleSaddle_Prob(out.uni1DSPASCC_FAST$root,X[G12,],muhat_m,u1,log.p=FALSE,method = "DSPASCC",pval.comp = "both",Var_U_beta = Var_U_beta,Var_U_star = Var_U_star,speedup = 1)
                p2DSPASCC_FAST1<-Get_DoubleSaddle_Prob(out.uni2DSPASCC_FAST$root,X[G12,],muhat_m,u1_inv,log.p=FALSE,method = "DSPASCC",pval.comp = "both",Var_U_beta = Var_U_beta,Var_U_star = Var_U_star,speedup = 1)

		p1DSPASCC_FAST2<-Get_DoubleSaddle_Prob(out.uni1DSPASCC_FAST$root,X[G12,],muhat_m,u1,log.p=FALSE,method = "DSPASCC",pval.comp = "both",Var_U_beta = Var_U_beta,Var_U_star = Var_U_star,speedup = 2)
                p2DSPASCC_FAST2<-Get_DoubleSaddle_Prob(out.uni2DSPASCC_FAST$root,X[G12,],muhat_m,u1_inv,log.p=FALSE,method = "DSPASCC",pval.comp = "both",Var_U_beta = Var_U_beta,Var_U_star = Var_U_star,speedup = 2)

		#pvalDSPA = p1DSPA + p2DSPA
		#pvalDSPAFCC = p1DSPAFCC + p2DSPAFCC
		pvalDSPASCC = p1DSPASCC + p2DSPASCC
		pvalSPA = p1SPA + p2SPA
		#pvalSPAFCC = p1SPAFCC + p2SPAFCC
		pvalSPASCC = p1SPASCC + p2SPASCC
		pvalDSPASCC_FAST1 = p1DSPASCC_FAST1 + p2DSPASCC_FAST1
		pvalDSPASCC_FAST2 = p1DSPASCC_FAST2 + p2DSPASCC_FAST2
		pvals[AddToRow,] = c(pvalSPA,pvalSPASCC,pvalDSPASCC,pvalDSPASCC_FAST1,pvalDSPASCC_FAST2)
		AddToRow = AddToRow + 1
 	}#end if Cutoff	

  }else{
  	SimsDone = sims - 1

  }


}#end for

pvals = pvals[1:(AddToRow-1),]
l = list(pvals,SimsDone)
save(l, file = paste("/home/pj268/project/paper1/TypeISims/","EUBRun",arID,"MAF",maf,"ccratio",ccratio,".RData",sep = ""))



