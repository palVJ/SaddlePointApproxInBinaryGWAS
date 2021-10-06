# 
# alpha = 0.05
# 
# n = 1000
# 
# n2 = 0
# n1 = 20
# n0 = 980
# 
# g = c(rep(0,n0),rep(1,n1),rep(2,n2))
# gtilde = g - mean(g)
# x = cbind(g,rep(1,n))
# 
# mu = 0.32 #two-sided SPA ok
# mu = 0.24 #two-sided SPA fails
# mu = 0.15 #one-sided SPA ok
# mu = 0.12 #one-sided SPA fails
# 
# mu = 0.88
# 
# ran = c(-sum(g)*mu, (1-mu)*sum(g))
# 
# uvals = seq(ran[1],ran[2],1)
# pm = pmf_condscore(uvals,n0,n1,n2,mu*n,mu)
# 


find_t1error = function(ran,uvals,pm,alpha){
    
    if(abs(ran[1])<abs(ran[2])){
        # Right tail wider than left tail
        
        tmp = uvals + ran[1]
        uid = which.max(tmp[tmp<=0]) #The largest positive value that gives a two-tailed test
        
        pval = pm[1] + sum(pm[(uid):length(pm)])
        
        if(pval > alpha){
            # Only observed score that yields a one-tailed test
            # can give rejection of the null hypothesis
            utailid = uid
            tmp = pval
            while(tmp > alpha){
                utailid = utailid + 1
                tmp = sum(pm[utailid:length(pm)])
            }
            
            crit_pos = utailid
            
        }else{
            # Observed score that yields two-tailed test
            # gives rejection of null hypothesis
            # Look first at observations on the right tail
            uposid = uid
            unegid = 1
            tmp = pval
            while(tmp < alpha){
                uposid = uposid - 1
                unegid = unegid + 1
                tmp = sum(pm[uposid:length(pm)]) + sum(pm[1:unegid])
            }
            uposid = uposid + 1
            unegid = unegid - 1
            tmp = sum(pm[uposid:length(pm)]) + sum(pm[1:unegid])
            
            crit_pos = uposid
            
        }
        
        # Look at observations at the left tail
        uneg = ran[1]
        upos = uneg-sign(uneg)*ceiling(2*abs(uneg))
        uposid = which(round(uvals,3) == round(upos,3))
        unegid = 1
        
        pval = pm[1] + sum(pm[uposid:length(pm)])
        if(pval<alpha){
            # Two-tailed test from an observation on the left tail can give rejection
            tmp = pval
            while(tmp < alpha){
                uposid = uposid-1
                unegid = unegid+1
                tmp = sum(pm[uposid:length(pm)]) + sum(pm[1:unegid])
            }
            uposid = uposid + 1
            unegid = unegid - 1
            tmp = sum(pm[uposid:length(pm)]) + sum(pm[1:unegid])
            
            crit_neg = unegid
            
            #The rejection region is
            #uvals[1:crit_neg] and
            #uvals[crit_pos:length(uvals)]
            
            t1error = sum(pm[1:crit_neg]) + 
                sum(pm[crit_pos:length(pm)])
            
        }else{
            crit_neg = NA
            
            #Rejection region is only 
            #uvals[crit_pos:length(uvals)]
            
            t1error = sum(pm[crit_pos:length(pm)])
        }
        
        return(c(t1error,crit_pos, crit_neg))
        
        
        
        
    }else if (abs(ran[1])>abs(ran[2])){
        # Left tail wider than right tail
        
        tmp = uvals + ran[2]
        uid = which.max(tmp[tmp<=0])+1 # den minste negative u-verdien som gir to-sidig test
        
        len = length(pm)
        
        pval = pm[len] + sum(pm[1:uid])
        
        if(pval > alpha){
            # Only negative score observation that yields one-tailed test
            # can give rejection
            utailid = uid
            tmp = pval
            while(tmp > alpha){
                utailid = utailid - 1
                tmp = sum(pm[1:utailid])
            }
            
            crit_neg = utailid
            
        }else{
            # Negative score observation that gives two-tailed test
            # can give rejection
            unegid = uid
            uposid = len
            tmp = pval
            while(tmp < alpha){
                uposid = uposid - 1
                unegid = unegid + 1
                tmp = sum(pm[uposid:length(pm)]) + sum(pm[1:unegid])
            }
            uposid = uposid + 1
            unegid = unegid - 1
            tmp = sum(pm[uposid:length(pm)]) + sum(pm[1:unegid])
            
            crit_neg = unegid
            
        }
        
        # Now we consider the right tail
        upos = ran[2]
        uneg = upos-sign(upos)*ceiling(2*abs(upos))
        unegid = which(round(uvals,3) == round(uneg,3))
        uposid = len
        
        pval = pm[len] + sum(pm[1:unegid])
        if(pval<alpha){
            # Two-tailed test can give rejection
            tmp = pval
            while(tmp < alpha){
                uposid = uposid-1
                unegid = unegid+1
                tmp = sum(pm[uposid:length(pm)]) + sum(pm[1:unegid])
            }
            uposid = uposid + 1
            unegid = unegid - 1
            tmp = sum(pm[uposid:length(pm)]) + sum(pm[1:unegid])
            
            crit_pos = uposid
            
            #Rejection region is
            #uvals[1:crit_neg] and
            #uvals[crit_pos:length(uvals)]
            
            t1error = sum(pm[1:crit_neg]) + 
                sum(pm[crit_pos:length(pm)])
            
        }else{
            crit_pos = NA
            
            #Rejection region is only
            #uvals[1:crit_neg]
            
            t1error = sum(pm[1:crit_neg])
        }
        
        return(c(t1error,crit_pos, crit_neg))
        
        
        
        
    }else{
        
        # Left and right tail is equally wide
        # We can only get two-tailed tests
        
        len = length(pm)
        pval = pm[1] + pm[len]
        if(pval > alpha){
            return(c(0,NA,NA))
        }else{
            tmp = pval
            unegid = 1
            uposid = len
            while(tmp < alpha){
                unegid = unegid + 1
                uposid = uposid - 1
                tmp = sum(pm[1:unegid]) + sum(pm[uposid:len])
            }
            unegid = unegid - 1
            uposid = uposid + 1
            tmp = sum(pm[1:unegid]) + sum(pm[uposid:len])
            
            return(c(tmp,uposid,unegid))
            
        }
        
        
    }
    
    
    
}


