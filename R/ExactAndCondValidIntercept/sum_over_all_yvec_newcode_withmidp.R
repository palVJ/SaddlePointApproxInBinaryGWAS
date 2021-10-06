

start = Sys.time()

library(ggplot2)
library(reshape2)
library(ggpubr)
library(viridis)

source("~/SaddlePointApproxInBinaryGWAS/R/ExactAndCondValidIntercept/CGFs.R")
source("~/SaddlePointApproxInBinaryGWAS/R/ExactAndCondValidIntercept/singlesadle.R")
source("~/SaddlePointApproxInBinaryGWAS/R/ExactAndCondValidIntercept/doublesadle.R")
source("~/SaddlePointApproxInBinaryGWAS/R/ExactAndCondValidIntercept/twosided_test.R")
source("~/SaddlePointApproxInBinaryGWAS/R/ExactAndCondValidIntercept/exact_test_intercept.R")
source("~/SaddlePointApproxInBinaryGWAS/R/ExactAndCondValidIntercept/exact_t1error.R")


# Compute exact probability of Type I error for each method: Exact, ESPA, ESPACC, DSPACC and normal approximation


# Significance level alpha = 0.05 and alpha = 5e-05 is chosen.
# n (sample size)
# Choose composition of g-vector (n0, n1, n2)

n = 1000
n2 = 0; n1 = 20; n0 = 980
g = c(rep(0,n0),rep(1,n1),rep(2,n2))
gtilde = g - mean(g)
x = cbind(g,rep(1,n))

for(r in 1:2){

if(r == 1){
    alpha = 0.05    
}else{
    alpha = 5e-05
}

ncases = 1:(n-1)


# In these vectors, the exact probability of type I error for each muhat
# and for each method is stored.

signExact = c()
signSPA = c()
signSPAcc = c()
signDSPAcc = c()
signNorm = c()

c = 0
#Evalute for muhat = 0.001 to 0.999
for(nc in ncases){
    c = c + 1
    
    print(nc)
    
    mu = nc/n    #muhat for this number of cases
    
    muhat = rep(mu,n) #muhat a constant for intercept model
    
    ran = round(c(-sum(g)*mu, (1-mu)*sum(g)),2) # The support can be no larger than this range
    uvals = round(seq(ran[1],ran[2],1),2) # alle possible score observations in this range
    
    pm = pmf_condscore(uvals,n0,n1,n2,nc,mu) # Exact probabilities for each observation
    #Correct the range as ran is lower and upper bounds, not exact bounds:
    #We can find the true range by looking at which uvals that have zero probability according to the exact distribution
    gtz = which(pm > 0)
    ran_true = c(uvals[gtz[1]],uvals[gtz[length(gtz)]])
    uvals = round(seq(ran_true[1],ran_true[2],1),2)
    
    # The function find_t1error finds the proability of type I error
    # for the exact test, in addition to upper t1[2] and lower t1[3]
    # critical value (NA if there can be no rejection)
    
    t1 = find_t1error(ran,uvals,pm,alpha)  
    
    signExact[c] = t1[1] # exact test type I error
    crit = t1[2] # Critical value in right tail from exact test
    
    
    ############################################################    
    if( !is.na(crit) & (((uvals[crit]-3) > abs(ran[1])) | (uvals[crit]-3) < 0)){
        # If TRUE:
        # The distribution to the test statistic is right-skew,
        # and we start the search three grid-values closer to zero from the critical point of exact test (rejection may occur before critical point for approximation methods),
        # for the approximation methods. Even at this start point, a one-tailed test will be calculated,
	#and we assume only a one-tailed test will give rejection
	#Also if starting point is negative, critical point is not far so zero, indicating a highly right-skew distribution, and we will assume only a one-tailed test will give rejection.
    ############################################################ 
        
        
        #################
        # ESPA
        #################
        uid = max(crit-3,which(uvals>0)[1]) 
        
        # uid is where we start our search, tre grid-values closer
        # to zero than exact test (or smallest positive u-value)
        
        # Calculate probability of rejection for uvals[uid]
        # and iterate along the tail until we find 
        # the critical value for ESPA
        
        tmpSPA = singlesaddle2_BN(uvals[uid],gtilde,muhat,ran)
        while((tmpSPA > alpha) & (uid < length(uvals_true))){
            uid = uid+1
            if(uid == length(uvals_true)){
                # NB: We need to do correction is we are on the boundary on the support (saddlepoint approximation not defined on the boundary of the support)
                tmpSPA = singlesaddle2_BN(uvals[uid]-0.5,gtilde,muhat,ran)
            }else{
                tmpSPA = singlesaddle2_BN(uvals[uid],gtilde,muhat,ran)
            }
        }
        if(tmpSPA > alpha){
            # Even on the boundary there was no rejection, probability of Type I error is zero
            signSPA[c] = 0
        }else{
            signSPA[c] = sum(pm[uid:length(pm)])
        }
        critSPA = uid
        
        
        #################
        # ESPACC
        #################
        uid = max(crit-3,which(uvals>0)[1])
        
        tmpSPAcc = singlesaddle2_SCC_BN(uvals[uid],gtilde,muhat,ran)
        while((tmpSPAcc > alpha) & (uid < length(uvals_true))){
            uid = uid+1
            if(uid == length(uvals_true)){
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[uid]-0.5,gtilde,muhat,ran)
            }else{
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[uid],gtilde,muhat,ran)
            }
        }
        
        if(tmpSPAcc > alpha){
            signSPAcc[c] = 0
        }else{
            signSPAcc[c] = sum(pm[uid:length(pm)])
        }
        critSPAcc = uid
        
        
        #################
        # DSPACC
        #################
        uid = max(crit-3,which(uvals>0)[1])
        
        tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uid],x,muhat,ran)
        while((tmpDSPAcc > alpha) & (uid < length(uvals_true))){
            uid = uid + 1
            if(uid == length(uvals_true)){
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uid]-0.5,x,muhat,ran)
            }else{
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uid],x,muhat,ran)
            }
        }
        
        if(tmpDSPAcc > alpha){
            signDSPAcc[c] = 0
        }else{
            signDSPAcc[c] = sum(pm[uid:length(pm)])
        }
        critDSPAcc = uid
        
        
    ############################################################    
    }else if(!is.na(crit)){
        # If TRUE:
        # Exact test has (at least) rejection at some place on the right tail
        # but there may be rejection on the left tail as well since the distribution is not that skewed
        # We need to investigate both tails:
    ############################################################    
        
        #################
        # ESPA
        #################
        
        # Explore the right tail first:
        # uposid is where we start the search
        uposid = max(crit-3,which(uvals>0)[1])
        
        tmpSPA = singlesaddle2_BN(uvals[uposid],gtilde,muhat,ran)
        while(tmpSPA > alpha & uposid < length(uvals)){
            uposid = uposid+1
            if(uposid == length(uvals)){
                tmpSPA = singlesaddle2_BN(uvals[uposid]-0.5,gtilde,muhat,ran)
            }else{
                tmpSPA = singlesaddle2_BN(uvals[uposid],gtilde,muhat,ran)
            }
        }
        
        if(tmpSPA > alpha){
            # No rejection at right tail
            critSPApos = NA
        }else{
            # Critical value for ESPA at right tail
            critSPApos = uposid
        }
        
        # Next, explore left tail
        # unegid is where we start the search
        # We start on the end of the tail and iterate
        # toward zero
        unegid = 1
        tmpSPA = singlesaddle2_BN(uvals[unegid]+0.5,gtilde,muhat,ran)
        if(tmpSPA < alpha){ 
        # If rejection at the end of the tail, continue to search towards zero
            while(tmpSPA < alpha){
                unegid = unegid+1
                tmpSPA = singlesaddle2_BN(uvals[unegid],gtilde,muhat,ran)
            }
            unegid = unegid - 1
            critSPAneg = unegid
            
            if(is.na(critSPApos)){
                # No rejection at right tail
                # Probability for type I error given by:
                signSPA[c] = sum(pm[1:critSPAneg])
            }else{
                # Rejections at both tails:
                signSPA[c] = sum(pm[1:critSPAneg]) + sum(pm[critSPApos:length(pm)])
            }
            
        }else{
        # Not possible to reject on the left tail.
            
            critSPAneg = NA
            if(is.na(critSPApos)){
                # No rejection, type I error probability is zero
                signSPA[c] = 0
            }else{
                # Type I error probability contribution from the right tail
                signSPA[c] = sum(pm[critSPApos:length(pm)])
            }
        }
        

        #################
        # ESPACC
        #################
        
        # From right tail
        uposid = max(crit-3,which(uvals>0)[1])
        tmpSPAcc = singlesaddle2_SCC_BN(uvals[uposid],gtilde,muhat,ran)
        
        while(tmpSPAcc > alpha & uposid < length(uvals)){
            uposid = uposid+1
            if(uposid == length(uvals)){
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[uposid]-0.5,gtilde,muhat,ran)
            }else{
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[uposid],gtilde,muhat,ran)
            }
        }
        
        if(tmpSPAcc > alpha){
            critSPAccpos = NA
        }else{
            critSPAccpos = uposid
        }
        
        # From left tail
        unegid = 1
        tmpSPAcc = singlesaddle2_SCC_BN(uvals[unegid]+0.5,gtilde,muhat,ran)
        if(tmpSPAcc < alpha){
            while(tmpSPAcc < alpha){
                unegid = unegid+1
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[unegid],gtilde,muhat,ran)
            }
            unegid = unegid-1
            critSPAccneg = unegid
            
            if(is.na(critSPAccpos)){
                signSPAcc[c] = sum(pm[1:critSPAccneg])
            }else{
                signSPAcc[c] = sum(pm[1:critSPAccneg]) + sum(pm[critSPAccpos:length(pm)])
            }
            
        }else{
            critSPAccneg = NA
            if(is.na(critSPAccpos)){
                signSPAcc[c] = 0
            }else{
                signSPAcc[c] = sum(pm[critSPAccpos:length(pm)])
            }
        }
        
        
        #################
        # DSPACC
        #################
        
        # From right tail
        uposid = max(crit-3,which(uvals>0)[1])
        tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uposid],x,muhat,ran)
        
        while(tmpDSPAcc > alpha & uposid < length(uvals)){
            uposid = uposid+1
            if(uposid == length(uvals)){
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uposid]-0.5,x,muhat,ran)
            }else{
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uposid],x,muhat,ran)
            }
        }
        if(tmpDSPAcc > alpha){
            critDSPAccpos = NA
        }else{
            critDSPAccpos = uposid
        }
        
        # From left tail
        unegid = 1
        tmpDSPAcc = doublesaddle2_SCC_BN(uvals[unegid]+0.5,x,muhat,ran)
        if(tmpDSPAcc < alpha){
            while(tmpDSPAcc < alpha){
                unegid = unegid+1
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[unegid],x,muhat,ran)
            }
            unegid = unegid-1
            critDSPAccneg = unegid
            
            if(is.na(critDSPAccpos)){
                signDSPAcc[c] = sum(pm[1:critDSPAccneg])
            }else{
                signDSPAcc[c] = sum(pm[1:critDSPAccneg]) + sum(pm[critDSPAccpos:length(pm)])
            }
            
        }else{
            critDSPAccneg = NA
            if(is.na(critDSPAccpos)){
                signDSPAcc[c] = 0
            }else{
                signDSPAcc[c] = sum(pm[critDSPAccpos:length(pm)])
            }
        }
        
        
    ############################################################    
    }else{
        # Distribution is left-skewed and rejection only in left
        # rail
    ############################################################ 
        crit = t1[3] # exact test critical value in left tail
        
        
        #################
        # ESPA
        #################
        
        # Start the seach in uid tre grid-value closer to zero from critical value
        # or least negative u-value
        # Seach towards the end of the left tail
        uid = min(crit+3,which(uvals<0)[length(which(uvals<0))])
        
        tmpSPA = singlesaddle2_BN(uvals[uid],gtilde,muhat,ran)
        while((tmpSPA > alpha) & (uid > 1)){
            uid = uid-1
            if(uid == 1){
                tmpSPA = singlesaddle2_BN(uvals[uid]+0.5,gtilde,muhat,ran)
            }else{
                tmpSPA = singlesaddle2_BN(uvals[uid],gtilde,muhat,ran)
            }
        }
        
        if (tmpSPA > alpha){
            # Found no rejection at left tail for ESPA 
            signSPA[c] = 0
        }else{
            signSPA[c] = sum(pm[1:uid])
        }
        critSPA = uid
        
        
        #################
        # SPACC
        #################
        uid = min(crit+3,which(uvals<0)[length(which(uvals<0))])
        
        tmpSPAcc = singlesaddle2_SCC_BN(uvals[uid],gtilde,muhat,ran)
        while(tmpSPAcc > alpha & uid > 1){
            uid = uid-1
            if(uid == 1){
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[uid]+0.5,gtilde,muhat,ran)
            }else{
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[uid],gtilde,muhat,ran)
            }
        }
        
        if (tmpSPAcc > alpha){
            signSPAcc[c] = 0
        }else{
            signSPAcc[c] = sum(pm[1:uid])
        }
        critSPAcc = uid
        
        
        #################
        # DSPACC
        #################
        uid = min(crit+3,which(uvals<0)[length(which(uvals<0))])
        
        tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uid],x,muhat,ran)
        while(tmpDSPAcc > alpha & uid > 1){
            uid = uid - 1
            if(uid == 1){
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uid]+0.5,x,muhat,ran)
            }else{
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uid],x,muhat,ran)
            }
        }
        
        if (tmpDSPAcc > alpha){
            signSSPAcc[c] = 0
        }else{
            signDSPAcc[c] = sum(pm[1:uid])
        }
        critDSPAcc = uid
        
        
        
    }
    
    
    ############################################################ 
    # Normal approximation
    ############################################################ 
    sdu = sqrt(mu*(1-mu)*sum(gtilde^2)) # standard deviation
    steps = round(min(uvals[uvals>0]),2) # The decimals of the positive grid points
    steps_neg = abs(steps-1) # The decimals of the negative integer grid points
    tmp = qnorm(0.5*alpha,0,sdu,lower.tail = FALSE) # (positive) critical value
    
    # Find the grid points on the support of the exact distribution closest to the positive and negative critical value of the normal approximaton
    
    # Right tail
    if(0 %in% uvals){
        crit_norm = ceiling(tmp)
        crit_norm_inv = floor(-tmp)
    }else if((tmp - floor(tmp)) < steps){
        crit_norm = floor(tmp) + steps
        crit_norm_inv = ifelse(-tmp < ceiling(-tmp) -steps_neg,ceiling(-tmp) -steps_neg-1,ceiling(-tmp) -steps_neg)
    }else{
        crit_norm = ceiling(tmp) + steps
        crit_norm_inv = ifelse(-tmp < ceiling(-tmp) -steps_neg,ceiling(-tmp) -steps_neg-1,ceiling(-tmp) -steps_neg)
    }
    
    
    if(crit_norm_inv>=ran[1] & crit_norm <= ran[2]){
        indpos = which(round(uvals,2) == round(crit_norm,2))
        indneg = which(round(uvals,2) == round(crit_norm_inv,2))
        signNorm[c] = sum(pm[1:indneg]) + sum(pm[indpos:length(pm)])
    }else if (crit_norm_inv<ran[1] & crit_norm <= ran[2]){
        indpos = which(round(uvals,2) == round(crit_norm,2))
        signNorm[c] = sum(pm[indpos:length(pm)])
    }else if (crit_norm_inv>=ran[1] & crit_norm > ran[2]){
        indneg = which(round(uvals,2) == round(crit_norm_inv,2))
        signNorm[c] = sum(pm[1:indneg])
    }else{
        signNorm[c] = 0
    }
}


signAll = cbind(ncases,
                signExact, 
                signSPA, 
                signSPAcc, 
                signDSPAcc, 
                signNorm,
                signMidP)


nmin = signAll[1,1]
nmax = signAll[dim(signAll)[1],1]

mumin = 0.05
mumax = 0.95
mus = round(seq(mumin,mumax,0.01),2)

# Type I error for each method for different values of true mu values (under H0)
t1error = matrix(NA,ncol = 6, nrow = length(mus))

# Probability of getting a data set where
# test is not valid (conditionally invalid test)
probab_invalid = matrix(NA,ncol = 6, nrow = length(mus))

muno = 0

for(mu in mus){
    muno = muno + 1
    yn1 = nmin:nmax # antall kasus vi skal se på
    pn1 = dbinom(yn1,n,mu) # sannsynlighet for antall kasus
    t1error[muno,] = c(sum(signAll[,2]*pn1),
                       sum(signAll[,3]*pn1),
                       sum(signAll[,4]*pn1),
                       sum(signAll[,5]*pn1),
                       sum(signAll[,6]*pn1),
                       sum(signAll[,7]*pn1))
    
    probab_invalid[muno,] = c(
        sum((signAll[,2]>alpha)*pn1),
        sum((signAll[,3]>alpha)*pn1),
        sum((signAll[,4]>alpha)*pn1),
        sum((signAll[,5]>alpha)*pn1),
        sum((signAll[,6]>alpha)*pn1),
        sum((signAll[,7]>alpha)*pn1))
}

minProb = apply(probab_invalid,2,min)
minProbID = apply(probab_invalid,2,which.min)
maxProb = apply(probab_invalid,2,max)
maxProbID = apply(probab_invalid,2,which.max)

prop_invalid = c(sum((signAll[,2]>alpha))/dim(signAll)[1],
                 sum((signAll[,3]>alpha))/dim(signAll)[1], 
                 sum((signAll[,4]>alpha))/dim(signAll)[1],
                 sum((signAll[,5]>alpha))/dim(signAll)[1],
                 sum((signAll[,6]>alpha))/dim(signAll)[1],
                 sum((signAll[,7]>alpha))/dim(signAll)[1])


which(signAll[,4]>alpha)
signAll[c(401,589),1]



colnames(t1error) = c("Exact","SPA","SPACC","DSPACC","Norm","MidP")
colnames(probab_invalid) = c("Exact","SPA","SPACC","DSPACC","Norm","MidP")

t1error = as.data.frame(t1error)
probab_invalid = as.data.frame(probab_invalid)


df= data.frame(mu = rep(mus,5))
df$t1errors = c(t1error$Exact,t1error$SPA,t1error$SPACC,t1error$DSPACC,t1error$Norm)
df$Method = c(rep("Exact",nrow(t1error)),rep("ESPA",nrow(t1error)),rep("ESPA-CC",nrow(t1error)),rep("DSPA-CC",nrow(t1error)),rep("Norm",nrow(t1error)))


df2= data.frame(mu = rep(mus,5))
df2$probab_invalids = c(probab_invalid$Exact,probab_invalid$SPA,probab_invalid$SPACC,probab_invalid$DSPACC,probab_invalid$Norm)
df2$Method = c(rep("Exact",nrow(t1error)),rep("ESPA",nrow(t1error)),rep("ESPA-CC",nrow(t1error)),rep("DSPA-CC",nrow(t1error)),rep("Norm",nrow(t1error)))

if(r == 1){
p1 = ggplot(df, aes(x = mu, y = t1errors, color = Method)) +
    theme_classic()+
    geom_line() +
    geom_abline(intercept = alpha,slope = 0) +
    scale_y_continuous("Overall type I error probabilities") +
    coord_cartesian(ylim=c(0, 0.063)) +
    scale_color_viridis(option="turbo",discrete = TRUE) +
    ggtitle(paste0("Signif. level ",alpha)) +
    labs(x=expression(mu))+
    scale_x_continuous(breaks = seq(0,1,0.1))

p1


p2 = ggplot(df2, aes(x = mu, y = probab_invalids,
                     color = Method)) +
    theme_classic()+
    geom_line() +
    coord_cartesian(ylim= c(0,1)) +
    scale_y_continuous("Probability of conditionally invalid test") +
    scale_color_viridis(option="turbo",discrete = TRUE) +
    ggtitle(paste0("Signif. level ",alpha)) +
    labs(x=expression(mu))+
    scale_x_continuous(breaks = seq(0.1,0.9,0.1)) #+ 


p2
}else{
    p3 = ggplot(df, aes(x = mu, y = t1errors, color = Method)) +
        theme_classic()+
        geom_line() +
        geom_abline(intercept = alpha,slope = 0) +
        scale_y_continuous("Overall type I error probabilities") +
        coord_cartesian(ylim=c(0, 1e-4)) +
        scale_color_viridis(option="turbo",discrete = TRUE) +
        ggtitle(paste0("Signif. level ",alpha)) +
        labs(x=expression(mu))+
        scale_x_continuous(breaks = seq(0,1,0.1))
    
    p3
    
    
    p4 = ggplot(df2, aes(x = mu, y = probab_invalids,
                         color = Method)) +
        theme_classic()+
        geom_line() +
        coord_cartesian(ylim= c(0,1)) +
        scale_y_continuous("Probability of conditionally invalid test") +
        scale_color_viridis(option="turbo",discrete = TRUE) +
        ggtitle(paste0("Signif. level ",alpha)) +
        labs(x=expression(mu))+
        scale_x_continuous(breaks = seq(0.1,0.9,0.1)) #+ 
    
    p4    
    
}

}

ggarrange(p1,p3,p2,p4,byrow=FALSE,nrow = 2,ncol = 2,common.legend = TRUE)



# plot(x = signAll[,1], y = signAll[,3], type = "l")
# lines(x = signAll[,1], y = signAll[,6], type = "l", col = "red")
# abline(h = alpha)
# 
# plot(x = signAll[,1], y = signAll[,6], type = "l")
# lines(x = signAll[,1], y = signAll[,3], type = "l", col = "red")
# abline(h = alpha)
# 
# 
# 
# minProb = c()
# maxProb = c()





