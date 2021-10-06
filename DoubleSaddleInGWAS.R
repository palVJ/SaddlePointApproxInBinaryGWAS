ScoreTest_NullModel_Get_X1 = function(X1)
{
        q1<-ncol(X1)
        if(q1>=2)
        {
                if(sum(abs(X1[,1]-X1[,2]))==0)
                {
                        X1=X1[,-2]
                        q1<-q1-1
                }
        }
        qr1<-qr(X1)
        if(qr1$rank < q1){

                X1.svd<-svd(X1)
                X1 = X1.svd$u[,1:qr1$rank]
        }

        return(X1)
}


ScoreTest_NULL_Model <- function(formula, data=NULL)
{
        X1<-model.matrix(formula,data=data)
        X1<-ScoreTest_NullModel_Get_X1(X1)

        glmfit= glm(formula, data=data, family = "binomial")
        convflag=0
        if(glmfit$converged)
        {
          muhat = glmfit$fitted.values
          if(mean(muhat)/mean(glmfit$y)>0.001 & (1-mean(muhat))/(1-mean(glmfit$y))>0.001)       convflag<-1     #Check that the null model converged properly with glm
        }
        if(convflag==0)
        {
          firthfit=fast.logistf.fit(x=X1,y=glmfit$y)
          eta<-X1%*%firthfit$beta
          muhat = as.vector(exp(eta)/(1+exp(eta)))
        }

        W = muhat*(1-muhat)
	      XW = t(X1 * W)
	      XWX = t(W*X1)%*%X1
        XWX_inv= solve(t(X1)%*%(X1 * W))
        XXVX_inv= X1 %*% XWX_inv

        re<-list(y=glmfit$y, muhat=muhat, W=W, X2=cov, XW=XW,XWX = XWX, XXWX_inv =XXWX_inv, convflag = convflag)
        class(re)<-"SA_NULL"
        return(re)
}

ScoreTest_NULL_Model_fast = function(pheno,cov){

	fastglmfit = fastglm::fastglm(y = pheno, x = cov,family = "binomial")
	convflag=0
  if(fastglmfit$converged){
    muhat = fastglmfit$fitted.values
    if(mean(muhat)/mean(fastglmfit$y)>0.001 & (1-mean(muhat))/(1-mean(fastglmfit$y))>0.001)       convflag<-1     #Check that the null model converged properly with glm

  }else{
	  stop("fastglm did not converge")
	}

	W = muhat*(1-muhat)
	XW = t(cov * W)
	XWX = t(W*cov)%*%cov
	XWX_inv = solve(XWX)
	XXWX_inv = cov%*%XWX_inv

	re<-list(y=fastglmfit$y, muhat=muhat, W=W, X2=cov, XW=XW,XWX = XWX, XXWX_inv =XXWX_inv, convflag = convflag)
	class(re)<-"SA_NULL"
	return(re)
}

KU=function(t,X,mu){
  
  X = as.matrix(X)
  Xmu<-mu%*%X
  etx<-exp(X%*%t)
  sum(log(1-mu*(1-etx)))-sum(t*as.vector(Xmu))
}

KUd1=function(t,X,mu){
  
  X = as.matrix(X)
  emtx<-exp(-X%*%t)
  apply(as.numeric(mu*(1/((1-mu)*emtx+mu)-1))*X,2,sum)
}

KUd2=function(t,X,mu)
{
  X = as.matrix(X)
  emtx<-exp(-X%*%t)
  t(as.numeric((mu*(1-mu)*emtx/((1-mu)*emtx+mu)^2))*X)%*%X
}

KU_approx=function(t,X_m,mu_m,Var_U_star)
{
  X_m = as.matrix(X_m)
  Xmu<-mu_m%*%X_m
  etx<-exp(X_m%*%t)
  sum(log(1-mu_m*(1-etx)))-sum(t*as.vector(Xmu)) + 0.5*t(t)%*%(Var_U_star)%*%t
}

KUd1_approx=function(t,X_m,mu_m,Var_U_star)
{
  X_m = as.matrix(X_m)
  emtx<-exp(-X_m%*%t)
  apply(as.numeric(mu_m*(1/((1-mu_m)*emtx+mu_m)-1))*X_m,2,sum) + Var_U_star%*%t
}

KUd2_approx=function(t,X_m,mu_m,Var_U_star)
{
  X_m = as.matrix(X_m)
  emtx<-exp(-X_m%*%t)
  t(as.numeric((mu_m*(1-mu_m)*emtx/((1-mu_m)*emtx+mu_m)^2))*X_m)%*%X_m + Var_U_star 
}


getrootDoubleSaddle_K1<-function(init,muhat,u1,u1_max,X,maxiter=300,reltol = 1e-8,method,Var_U_star = NULL){


#Double SPA assuming continuous distribution or double spa with first continuity correction
if(method == "DSPA" | method == "DSPAFCC"){

	#check if u1 is equal to maximum value. If so, do continuity correction here:
	if(u1 == u1_max){
		u1 = u1 - 0.5
	}

	if(method == "DSPAFCC" && u1 < 0){
	#We compute p-value from the left tail meaning we want Pr(U1 =< u1) = 1 - Pr(U >= u1+1). So by default compute saddlepoint at u1+1 
		u1 = u1 + 1
	#Find saddlepoint either exactly or approximately but faster:
		if(is.null(Var_U_star)){
			optimize = optim(init,method = "BFGS",control = list(maxit = maxiter,reltol = reltol),function(t)KU(t,X,muhat)-u1*t[1], gr = function(t)KUd1(t,X,muhat)- c(u1,rep(0,length(init)-1)))
		}else{
			optimize = optim(init,method = "BFGS",control = list(maxit = maxiter,reltol = reltol),function(t)KU_approx(t,X,muhat,Var_U_star)-u1*t[1], gr = function(t)KUd1_approx(t,X,muhat,Var_U_star)- c(u1,rep(0,length(init)-1)))
		}
		return(list(root=optimize$par,convergence = optimize$convergence))

	}else{
	#We use DSPA or DSPAFCC with p-value computation at right tail.
		if(is.null(Var_U_star)){
			optimize = optim(init,method = "BFGS",control = list(maxit = maxiter,reltol = reltol),function(t)KU(t,X,muhat)-u1*t[1], gr = function(t)KUd1(t,X,muhat)- c(u1,rep(0,length(init)-1)))
		}else{
			optimize = optim(init,method = "BFGS",control = list(maxit = maxiter,reltol = reltol),function(t)KU_approx(t,X,muhat,Var_U_star)-u1*t[1], gr = function(t)KUd1_approx(t,X,muhat,Var_U_star)- c(u1,rep(0,length(init)-1)))
		}
        	return(list(root=optimize$par,convergence = optimize$convergence))
	}

}else{
#Second continuity correction meaning method == "DSPASCC"
	if(u1 > 0){
	#Compute offset saddlepoint:
		u1 = u1-0.5
		if(is.null(Var_U_star)){
			optimize = optim(init,method = "BFGS",control = list(maxit = maxiter,reltol = reltol),function(t)KU(t,X,muhat)-u1*t[1], gr = function(t)KUd1(t,X,muhat)- c(u1,rep(0,length(init)-1)))
		}else{
			optimize = optim(init,method = "BFGS",control = list(maxit = maxiter,reltol = reltol),function(t)KU_approx(t,X,muhat,Var_U_star)-u1*t[1], gr = function(t)KUd1_approx(t,X,muhat,Var_U_star)- c(u1,rep(0,length(init)-1)))
		}
        	return(list(root=optimize$par,convergence = optimize$convergence))	
	}

	else{
	#We compute p-value from the left tail meaning we want Pr(U =< u1) = 1 - Pr(U >= u1+1). So by default compute saddlepoint at u1+1
		u1 = u1+ 0.5
		if(is.null(Var_U_star)){
                	optimize = optim(init,method = "BFGS",control = list(maxit = maxiter,reltol=reltol),function(t)KU(t,X,muhat)-u1*t[1], gr = function(t)KUd1(t,X,muhat)- c(u1,rep(0,length(init)-1)))
		}else{
			optimize = optim(init,method = "BFGS",control = list(maxit = maxiter,reltol = reltol),function(t)KU_approx(t,X,muhat,Var_U_star)-u1*t[1], gr = function(t)KUd1_approx(t,X,muhat,Var_U_star)- c(u1,rep(0,length(init)-1)))
		}
                return(list(root=optimize$par,convergence = optimize$convergence))	

	}


}

}


getrootSaddleContCorr_K1<-function(init,muhat,q,u1,u1_max,G_tilde,maxiter=300,reltol = 1e-8,method,Var_U_star = NULL){

#First continuity correction
if(method == "SPAFCC"){

	#check if u1 is equal to maximum value. If so, do continuity correction here:
        if(u1 == u1_max){
                q = q - 0.5
        }

	if(u1 <0){
	#We compute p-value from the left tail meaning we want Pr(U =< u1) = 1 - Pr(U >= u1+1). So by default compute saddlepoint at u1+1, meaning q = q +1
                sp = SPAtest:::getroot_K1(init, mu=muhat, g=G_tilde, q=q+1)$root
		return(list(root = sp,convergence = 0))
	}else{
		sp = SPAtest:::getroot_K1(init, mu=muhat, g=G_tilde, q=q)$root
		return(list(root = sp,convergence = 0))
	}

}else{
#method is SPASCC

	if(u1< 0){
		sp = SPAtest:::getroot_K1(init, mu=muhat, g=G_tilde, q=q+0.5)$root
		return(list(root = sp,convergence = 0))
	}else{
		sp = SPAtest:::getroot_K1(init, mu=muhat, g=G_tilde, q=q-0.5)$root
                return(list(root = sp,convergence = 0))
	}

}



}

Get_Saddle_Prob_SPA<-function(that,muhat,G_tilde,u1,q,log.p=FALSE,pval.comp){

	k1<-SPAtest:::Korg(that, muhat, G_tilde)
	k2<-SPAtest:::K2(that, muhat, G_tilde)
	
	if(is.finite(k1) && is.finite(k2)){
		temp1 <- that*q - k1

	
		what<-sign(that)*(2*temp1)^{1/2}
		vhat<- that*(k2)^{1/2}
	
		if(pval.comp == "BN"){ 
			Z.test<-what + 1/what * log(vhat/what)	
	

			if(Z.test > 0){
				pval<-pnorm(Z.test, lower.tail = FALSE,log.p=log.p)
			}else{
				pval= pnorm(Z.test, lower.tail = TRUE,log.p=log.p)
			}
		}else if(pval.comp == "LR"){
		#Do Skovgaard p-value computation

                	S.test = pnorm(what)+dnorm(what)*((1/what)-(1/vhat))

                        if(u1 > 0){
                        	pval = 1 - S.test
                        }else{
                        	pval = S.test
                        }

		}else{
		#pval.comp == "both"

                	Z.test<-what + 1/what * log(vhat/what)


                        if(Z.test > 0){
                        	pvalBN = pnorm(Z.test, lower.tail = FALSE,log.p=log.p)
                        }else{
                        	pvalBN = pnorm(Z.test, lower.tail = TRUE,log.p=log.p)
                        }

                        S.test = pnorm(what)+dnorm(what)*((1/what)-(1/vhat))
                        if(u1 > 0){
                        	pvalLR = 1 - S.test
                        }else{
                        	pvalLR = S.test
                        }

                        pval = c(pvalLR,pvalBN)

              }
	
	}else{
		if(log.p){
			pval<- -Inf
		}else{
		pval<-0
		}
	}
	
	return(pval)
}


Get_DoubleSaddle_Prob<-function(that,X,muhat,u1,log.p=FALSE,method,pval.comp,Var_U_star = NULL,Var_U_beta = NULL,speedup = NULL){
        if(is.na(that[1])){
                return(0)
        }else{
		#compute CGFs. Either exact or approximately but faster:
		if(is.null(Var_U_star)){
                	k <-KU(that,X, muhat)
                	k2<-KUd2(that,X, muhat)
                	Det_Hess_that = det(k2)

			#The determinant of KUd2 under H_0: t(X2)%*%(pheno-muhat)=0 with saddlepoint that = 0. The determinant is with respect to X2 (not including G_tilde):
                	Det_Hess_t_tilde_hat = det(KUd2(rep(0,ncol(X)-1),X[,-1],muhat))
		}else{
		  k <-KU_approx(that,X, muhat,Var_U_star)
      k2<-KUd2_approx(that,X, muhat,Var_U_star)
      Det_Hess_that = det(k2)

			#The determinant of KUd2 under H_0: t(X2)%*%(pheno-muhat)=0 with saddlepoint that = 0. The determinant is with respect to X2 (not including G_tilde):
			#Dependent on which speedup:
      Det_Hess_t_tilde_hat = ifelse(speedup == 1,det(Var_U_beta),det(KUd2_approx(rep(0,ncol(X)-1),X[,-1], muhat,Var_U_star[2:ncol(X),2:ncol(X)])))
		}

                if(is.finite(k) && is.finite(k2)){

			if(method == "DSPA"){
				what = sign(that[1])*sqrt(2*(-k+that[1]*u1))
				vhat =  that[1]*sqrt(Det_Hess_that/Det_Hess_t_tilde_hat)

			}else if(method == "DSPAFCC"){

				if(u1 <0){
				#We want to compute p-value from the left tail meaning we want Pr(U1 =< u1) = 1 - Pr(U1 >= u1+1). So compute by default what and vhat at u1+1
                                        what = sign(that[1])*sqrt(2*(-k+that[1]*(u1+1)))										
				}else{
					what = sign(that[1])*sqrt(2*(-k+that[1]*u1))
				}
				vhat = (1-exp(-that[1]))*sqrt(Det_Hess_that/Det_Hess_t_tilde_hat)
			}else{
			#second continuity correction:

				if(u1 < 0){
                                #We want to compute p-value from the left tail meaning we want Pr(U1 =< u1) = 1 - Pr(U1 >= u1+1). So compute by default what at u1+1
                                        what = sign(that[1])*sqrt(2*(-k+that[1]*(u1+0.5)))
                                }else{
                                        what = sign(that[1])*sqrt(2*(-k+that[1]*(u1-0.5)))
                                }
				vhat = 2*sinh(that[1]/2)*sqrt(Det_Hess_that/Det_Hess_t_tilde_hat)
			}

			#Do Barndorff-Nielsen p-value computation
			if(pval.comp == "BN"){

                        	Z.test<-what + 1/what * log(vhat/what)


                        	if(Z.test > 0){
                                	pval = pnorm(Z.test, lower.tail = FALSE,log.p=log.p)
                        	}else{
                                	pval = pnorm(Z.test, lower.tail = TRUE,log.p=log.p)
                        	}


			}else if(pval.comp == "LR"){
			#Do Skovgaard p-value computation

				S.test = pnorm(what)+dnorm(what)*((1/what)-(1/vhat)) 

				if(u1 > 0){
					pval = 1 - S.test
				}else{
					pval = S.test
				}			

			}else{
			#pval.comp == "both"
			
				Z.test<-what + 1/what * log(vhat/what)


                                if(Z.test > 0){
                                        pvalBN = pnorm(Z.test, lower.tail = FALSE,log.p=log.p)
                                }else {
                                        pvalBN = pnorm(Z.test, lower.tail = TRUE,log.p=log.p)
                                }

				S.test = pnorm(what)+dnorm(what)*((1/what)-(1/vhat))
                                if(u1 > 0){
                                        pvalLR = 1 - S.test
                                }else{
                                        pvalLR = S.test
                                }
				
				pval = c(pvalLR,pvalBN)

			}

                }else{
                        if(log.p){
                                pval<- -Inf
                        }else{
                       		pval<-0
                        }
                     }

        	return(pval)
        }
}


Get_Saddle_ProbContCor<-function(that,G_tilde,muhat,u1,q,log.p=FALSE,method,pval.comp){

if(is.na(that[1])){
                return(0)
        }else{

		k<-SPAtest:::Korg(that, muhat, G_tilde)
		k2<-SPAtest:::K2(that,muhat,G_tilde)
                
		if(is.finite(k) && is.finite(k2)){


			if(method == "SPAFCC"){
				if(u1 <0){
                                #We want to compute p-value from the left tail meaning we want Pr(U1 <= u1) = 1 - Pr(U1 >= u1+1). So compute by default what and vhat at u1+1, or q = q + 1
                                        what = sign(that)*sqrt(2*(-k+that*(q+1)))
                                }else{
                                        what = sign(that)*sqrt(2*(-k+that*q))
                                }
 
                                vhat = (1-exp(-that))*(k2)^{1/2}
	
			}else{
			#second continuity correction
				if(u1 < 0){
                                #We want to compute p-value from the left tail meaning we want Pr(U1 <= u1) = 1 - Pr(U >= u1+1). So compute by default what at u1+1, or q = q + 1
                                        what = sign(that)*sqrt(2*(-k+that*(q+0.5)))
                                }else{
                                        what = sign(that)*sqrt(2*(-k+that*(q-0.5)))
                                }
                             
                                vhat = 2*sinh(that/2)*(k2)^{1/2}

			}

			#Do Barndorff-Nielsen p-value computation
                        if(pval.comp == "BN"){

                                Z.test<-what + 1/what * log(vhat/what)
				

                                if(Z.test > 0){
                                        pval = pnorm(Z.test, lower.tail = FALSE,log.p=log.p)
                                }else{
                                        pval = pnorm(Z.test, lower.tail = TRUE,log.p=log.p)
                                }


                        }else if(pval.comp == "LR"){
                        #Do Skovgaard p-value computation

                                S.test = pnorm(what)+dnorm(what)*((1/what)-(1/vhat))

                                if(u1 > 0){
                                        pval = 1 - S.test
                                }else{
                                        pval = S.test
                                }

                        }else{
                        #pval.comp == "both"

                                Z.test<-what + 1/what * log(vhat/what)


                                if(Z.test > 0){
                                        pvalBN = pnorm(Z.test, lower.tail = FALSE,log.p=log.p)
                                }else {
                                        pvalBN = pnorm(Z.test, lower.tail = TRUE,log.p=log.p)
                                }

                                S.test = pnorm(what)+dnorm(what)*((1/what)-(1/vhat))
                                if(u1 > 0){
                                        pvalLR = 1 - S.test
                                }else{
                                        pvalLR = S.test
                                }

                                pval = c(pvalLR,pvalBN)

                        }



		}else{
                        if(log.p){
                                pval<- -Inf
                        }else{
                                pval<-0
                        }
                     }

                return(pval)

	     }



}


DoubleSaddle_Prob<-function(u, u1_min, u1_max, muhat, X, Cutoff=2,alpha,output="P",nodes.fixed,nodes.init,log.p=FALSE,method,pval.comp)
{
        Ncovs = ncol(X)
        G_tilde = X[,1]

        var1<-sum(muhat * (1-muhat) * G_tilde^2)
        p1=NULL
        p2=NULL

        ObservedU1 = u[1]
        MinusObservedU1 = -u[1]

        # Noadj
        pval.noadj<-pchisq((ObservedScore)^2/var1, lower.tail = FALSE, df=1,log.p=log.p)
        convergence=TRUE


        if(Cutoff=="BE"){
                rho<-sum(((abs(g))^3)*muhat*(1-muhat)*(muhat^2+(1-muhat)^2))
                B<-0.56*rho*var1^(-3/2)
                p<-B+alpha/2
                Cutoff=ifelse(p>=0.496,0.01,qnorm(p,lower.tail=F))
        }

	if(Cutoff < 10^-1){
        	Cutoff=10^-1
        }


        if(abs(ObservedU1)/sqrt(var1) < Cutoff ){

                pval=pval.noadj
        }else{

		out.uni1 = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat,u,u1_max,X,method)

                if(MinusObservedU1 > u1_max || MinusObservedU1 < u1_min){
                #if minus observed score is outside range of density, no saddlepoint exists and we do one-sided p-value calculations
                        out.uni2 = list(root = NA,convergence = 0)
                }else{
                        out.uni2 = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat,-u,u1_max,X,method)
                }

		if(out.uni1$convergence==0 && out.uni2$convergence==0)
                {
                        p1<-tryCatch(Get_DoubleSaddle_Prob(out.uni1$root,X,muhat,u,log.p=FALSE,method,pval.comp),error=function(e) {
                                if(log.p) return(pval.noadj-log(2))
                                else return(pval.noadj/2)})
                        p2<-tryCatch(Get_DoubleSaddle_Prob(out.uni2$root,X,muhat,-u,log.p=FALSE,method,pval.comp),error=function(e) {
                                if(log.p) return(pval.noadj-log(2))
                                else return(pval.noadj/2)})
                        
			pval = p1+p2
                        convergence= 0
                }else {
                        print("Error_Converge")
                        pval<-pval.noadj
                        convergence= 2
                }
        }

	if(pval!=0 && pval.noadj/pval>10^3)
        {
                return(DoubleSaddle_Prob(u, u1_min, u1_max, muhat, X, Cutoff=Cutoff*2,alpha,output="P",nodes.fixed,nodes.init,log.p=FALSE,method,pval.comp))
        }else {
                return(list(p.value=pval, p.value.NA=pval.noadj,convergence=convergence, Score=ObservedU1/sqrt(var1)))
        }


}

TestDoubleSPA<-function(G, obj.null, Cutoff=2,alpha,output,nodes.fixed=NULL,nodes.init,log.p=FALSE,method,pval.comp)
{
        if(class(obj.null) != "SA_NULL"){
                stop("obj.null should be a returned object from ScoreTest_wSaddleApprox_NULL_Model")
        }

        pheno = obj.null$y
        muhat = obj.null$muhat

        n.g<-sum(G)
        if(n.g/(2*length(G))>0.5)
        {
                G<-2-G
                n.g<-sum(G)
        }

        G_tilde = G - (obj.null$XXVX_inv)%*%((obj.null$XV)%*%G)
        X = cbind(G_tilde,obj.null$X2)
        u = crossprod(as.matrix(X),as.matrix(pheno-muhat))
	u1_min = -sum(G*muhat)
        u1_max = sum(G*(1-muhat))

        out = DoubleSaddle_Prob(u,u1_min,u1_max,muhat, X = X, Cutoff=Cutoff,alpha,output="P",nodes.fixed,nodes.init,log.p=FALSE,method,pval.comp)

        return(out)
}

ScoreTest <-function(genos,pheno,cov,obj.null,minmac=5,Cutoff=2,alpha=5*10^-8,beta.out=FALSE,beta.Cutoff=5*10^-7,log.p=FALSE,pval.comp,methods){

        if(missing(obj.null)){

          if(missing(cov) || is.null(cov)){
            cov<-rep(1,length(pheno))
          }
          
          obj.null<-ScoreTest_NULL_Model_fast(pheno,cov)
        }
	
	
     
        if(is.vector(genos)){
		      m = 1
		      genos = as.matrix(genos)
	      }else{
		      m <- ncol(genos)		
	      }
     
        p.values<- matrix(data= NA, ncol = length(methods),nrow = m)
	      colnames(p.values) = methods
	      rownames(p.values) = colnames(genos)
	      Ninds = nrow(cov)
        for (i in 1:m){

		      G = genos[,i]
          ina<-which(is.na(G))

	
		      if(length(ina)>0){
            G[ina] <-median(G,na.rm=TRUE)
          }
          MAC<-min(sum(G),sum(2-G))
		
          if(MAC>=minmac){
                  
            muhat = obj.null$muhat
        		G_tilde = G - (obj.null$XXWX_inv)%*%((obj.null$XW)%*%G)
        		X = cbind(G,cov)
            Ncovs = ncol(X)
        		u = crossprod(as.matrix(X),as.matrix(pheno-muhat))
			      u1 = u[1]
        		var1<-sum(muhat*(1-muhat)*G_tilde^2)
        		score = u1/sqrt(var1)
        		u1_min = -sum(G*muhat)
        		u1_max = sum(G*(1-muhat))

			      pval_noraprx = pchisq(score^2, lower.tail = FALSE, df=1,log.p=log.p)

			      #If normal approximation "noraprx" is part of methods to use:
			      if("noraprx" %in% methods){
				      p.values[i,which(methods == "noraprx")] = pval_noraprx
			      }

			      #If score sufficiently large to do saddlepoint approximations:
			      if(abs(score)>= Cutoff){

              m1 = sum(muhat*G_tilde)
              q = sum(G_tilde*pheno)
			
				      #Find grid point closest to -u[1]
              u1_inv = u1 - sign(u1)*ceiling(abs(2*u1))
				      qinv = u1_inv + m1


				      if("SPA" %in% methods){
				
					      out.uni1SPA<-SPAtest:::getroot_K1(0, mu=muhat, g=G_tilde, q=q)
                out.uni2SPA<-SPAtest:::getroot_K1(0, mu=muhat, g=G_tilde, q=qinv)
				        #out.uni1SPA = optimize(function(t)KU(t,as.matrix(G_tilde),muhat)-u1*t,c(-10,10),tol=1e-10)$minimum
					      p1SPA = Get_Saddle_Prob_SPA(that = out.uni1SPA$root,u1 = u1, muhat = muhat, G_tilde = G_tilde, q = q,log.p=FALSE,pval.comp)
                p2SPA = Get_Saddle_Prob_SPA(that = out.uni2SPA$root,u1 = u1_inv, muhat = muhat, G_tilde = G_tilde, q = qinv,log.p=FALSE,pval.comp)				

					      p.values[i,which(methods == "SPA")] = p1SPA + p2SPA

				      }

              if("SPASCC" %in% methods){
					      out.uni1SPASCC = getrootSaddleContCorr_K1(0,muhat,q,u1,u1_max,G_tilde,method = "SPASCC")

					      if(-u1 > u1_max || -u1 < u1_min){
                #if minus observed u1 is outside range of density, no saddlepoint exists and we do one-sided p-value calculations
                  out.uni2SPASCC = list(root = NA,convergence = 0)
                }else{
            
                out.uni2SPASCC = getrootSaddleContCorr_K1(0,muhat,qinv,u1_inv,u1_max,G_tilde,method = "SPASCC")
                }

					      p1SPASCC = Get_Saddle_ProbContCor(out.uni1SPASCC$root,G_tilde,muhat,u1,q,log.p=FALSE,method = "SPASCC",pval.comp)
                p2SPASCC = Get_Saddle_ProbContCor(out.uni2SPASCC$root,G_tilde,muhat,u1_inv,qinv,log.p=FALSE,method = "SPASCC",pval.comp)
					      p.values[i,which(methods == "SPASCC")] = p1SPASCC + p2SPASCC

				      }

				      if("DSPASCC" %in% methods){

					      #Compute p-value DSPASCC:
					      out.uni1DSPASCC = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat,u1,u1_max,X,method = "DSPASCC")

					      if(-u1 > u1_max || -u1 < u1_min){
                #if minus observed U1 is outside range of density, no saddlepoint exists and we do one-sided p-value calculations
                  out.uni2DSPASCC = list(root = NA,convergence = 0)                                          
                }else{
                  out.uni2DSPASCC = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat,u1_inv,u1_max,X,method = "DSPASCC")
                }

					      p1DSPASCC<-Get_DoubleSaddle_Prob(out.uni1DSPASCC$root,X,muhat,u1,log.p=FALSE,method = "DSPASCC",pval.comp)
                p2DSPASCC<-Get_DoubleSaddle_Prob(out.uni2DSPASCC$root,X,muhat,u1_inv,log.p=FALSE,method = "DSPASCC",pval.comp)
					      p.values[i,which(methods == "DSPASCC")] = p1DSPASCC + p2DSPASCC

				      }

				      if(sum(grepl("DSPASCC_FAST",methods))>0){
					
				        #Compute variables for approximate computations of CGFs
                Var_U_beta = obj.null$XWX
                G12 = which(G > 0)
                muhat_m = muhat[G12]
                cov_m = cov[G12,]
                Var_U_star = matrix(0,ncol = Ncovs,nrow = Ncovs)
                Var_U_star[2:Ncovs,2:Ncovs] = Var_U_beta - t(muhat_m*(1-muhat_m)*cov_m)%*%cov_m
					      out.uni1DSPASCC_FAST = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat_m,u1,u1_max,X[G12,],method = "DSPASCC",Var_U_star = Var_U_star)

					    if(-u1 > u1_max || -u1 < u1_min){
                #if minus observed U1 is outside range of density, no saddlepoint exists and we do one-sided p-value calculations
                out.uni2DSPASCC_FAST = list(root = NA,convergence = 0)
              }else{
                out.uni2DSPASCC_FAST = getrootDoubleSaddle_K1(rep(0,Ncovs),muhat_m,u1_inv,u1_max,X[G12,],method = "DSPASCC",Var_U_star = Var_U_star)
              }
					
					    if("DSPASCC_FAST1" %in% methods){
					      p1DSPASCC_FAST1<-Get_DoubleSaddle_Prob(out.uni1DSPASCC_FAST$root,X[G12,],muhat_m,u1,log.p=FALSE,method = "DSPASCC",pval.comp,Var_U_beta = Var_U_beta,Var_U_star = Var_U_star,speedup = 1)
                p2DSPASCC_FAST1<-Get_DoubleSaddle_Prob(out.uni2DSPASCC_FAST$root,X[G12,],muhat_m,u1_inv,log.p=FALSE,method = "DSPASCC",pval.comp,Var_U_beta = Var_U_beta,Var_U_star = Var_U_star,speedup = 1)
						    p.values[i,which(methods == "DSPASCC_FAST1")] = p1DSPASCC_FAST1 + p2DSPASCC_FAST1
					    }


					    if("DSPASCC_FAST2" %in% methods){
					      p1DSPASCC_FAST2<-Get_DoubleSaddle_Prob(out.uni1DSPASCC_FAST$root,X[G12,],muhat_m,u1,log.p=FALSE,method = "DSPASCC",pval.comp,Var_U_beta = Var_U_beta,Var_U_star = Var_U_star,speedup = 2)
                p2DSPASCC_FAST2<-Get_DoubleSaddle_Prob(out.uni2DSPASCC_FAST$root,X[G12,],muhat_m,u1_inv,log.p=FALSE,method = "DSPASCC",pval.comp,Var_U_beta = Var_U_beta,Var_U_star = Var_U_star,speedup = 2)
						    p.values[i,which(methods == "DSPASCC_FAST2")] = p1DSPASCC_FAST2 +p2DSPASCC_FAST1
					    }
				    }
                                        

			     }else{
				    #Normal approximation is sufficient
				    p.values[i,] = pval_noraprx
			     }
			      
		      }

        }
	
        return(p.values)
}

GWAS_screen = function(genos,cov,pheno,obj.null,minmac=1,alpha=c(5e-5,5e-8),arID,outfilename){
	
	#The name of the snps:
  snp_names = colnames(genos)

	snps_of_interest = matrix(ncol = length(alpha),nrow = 0)
	

	if(is.vector(genos)){
                m = 1
                genos = as.matrix(genos)
        }else{
                m <- ncol(genos)
        }



	for(i in 1:ncol(genos)){

                G = unlist(genos[,..i])
                ina<-which(is.na(G))

                if(length(ina)>0){
                        G[ina] <-median(G,na.rm=TRUE)
                }
                MAC<-min(sum(G),sum(2-G))

                if(MAC>=minmac){

                        muhat = obj.null$muhat
                        G_tilde = G - (obj.null$XXWX_inv)%*%((obj.null$XW)%*%G)
                        X = cbind(G_tilde,cov)
                        Ncovs = ncol(X)
                        u = crossprod(as.matrix(X),as.matrix(pheno-muhat))
                        u1 = u[1]
                        var1<-sum(muhat*(1-muhat)*G_tilde^2)
                        score = u1/sqrt(var1)

                        pval_noraprx = pchisq(score^2, lower.tail = FALSE, df=1)

                        #If normal approximated p-value is less than alpha, collect name of snp
                        lessalpha = pval_noraprx < alpha
			if(sum(lessalpha) > 0){
				
				snp = c(NA,NA)
				snp[which(lessalpha)] = snp_names[i] 
				snps_of_interest = rbind(snps_of_interest,snp)

			}

		}
	}


	library(data.table)
	#write snp_names to file:
	for(j in 1:length(alpha)){

		snps = snps_of_interest[,j]
		snps = as.data.frame(snps[!is.na(snps)])
		fwrite(snps,paste(outfilename,alpha[j],"arID",arID,".txt",sep = ""),quote = FALSE,row.names = FALSE, col.names = FALSE,append = TRUE)
	}

	return(snps_of_interest)

}


