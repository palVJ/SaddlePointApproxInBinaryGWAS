

# Exact test, two-sided, intercept model

pmf_condscore = Vectorize(function(u,n0,n1,n2,np,mu){
    ustar = u + mu*(n1+2*n2)
    tmp = 0
    lognevner = lchoose(n,np)  
    for(k in max(ceiling(.5*(ustar-n1)),0):min(floor(.5*ustar),n2)){
        tmp = tmp + exp(lchoose(n0,np-ustar+k)+lchoose(n1,ustar-2*k)+lchoose(n2,k)-lognevner)
    }
    tmp
}, vectorize.args = "u")


exact_intercept_2 = function(u,n0,n1,n2,np,mu,ran){
    uvals = seq(ran[1],ran[2],1)
    pmfu = pmf_condscore(uvals,n0,n1,n2,np,mu)
    
    if(u < 0){
        upos = u-sign(u)*ceiling(2*abs(u))
        
        ind_pos = which(round(uvals,3) == round(upos,3))
        ind_neg = which(round(uvals,3) == round(u,3))
        
        pu = sum(pmfu[ind_pos:length(pmfu)])
        pl = sum(pmfu[1:ind_neg])
        
        return(pu + pl)
    }else if(u > abs(ran[1])){
        # Only right-tail
        ind_pos = which(round(uvals,3) == round(u,3))
        pu = sum(pmfu[ind_pos:length(pmfu)])
        return(pu)
    }else{
        uneg = u-sign(u)*ceiling(2*abs(u))
        
        ind_pos = which(round(uvals,3) == round(u,3))
        ind_neg = which(round(uvals,3) == round(uneg,3))
        
        pu = sum(pmfu[ind_pos:length(pmfu)])
        pl = sum(pmfu[1:ind_neg])
        
        return(pu + pl)
    }
}


