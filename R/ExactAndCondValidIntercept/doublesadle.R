
# Estimate P(Ugamma >= u | Ubeta = 0)

# u is the observed score (Ugamma)
# x is a covariate matrix with genotype-vector as the first column 
# and always vector of ones as the second column
# muhat is a vector of estimated means under H0


pp_doublesaddle_SCC_LR = function(u,x,muhat){
    x = as.matrix(x)
    nx = dim(x)[2]
    x2 = x[,2:nx]
    
    uoff = u - 0.5

    that = optim(rep(0,nx),function(t)k(t,x,muhat)-uoff*t[1])$par
    
    what = sign(that[1])*sqrt(2*(-k(that,x,muhat)+that[1]*uoff)) 
    
    vhat = 2*sinh(that[1]/2)*sqrt(det(D2k(that,x,muhat))/det(D2k(rep(0,nx-1),x2,muhat)))
    
    1-pnorm(what)-dnorm(what)*(1/what-1/vhat)
}



pp_doublesaddle_SCC_BN = function(u,x,muhat){
    x = as.matrix(x)
    nx = dim(x)[2]
    x2 = x[,2:nx]
    
    uoff = u - 0.5
    
    that = optim(rep(0,nx),function(t)k(t,x,muhat)-uoff*t[1])$par
    
    what = sign(that[1])*sqrt(2*(-k(that,x,muhat)+that[1]*uoff)) 
    
    vhat = 2*sinh(that[1]/2)*sqrt(det(D2k(that,x,muhat))/det(D2k(rep(0,nx-1),x2,muhat)))
    
    1 - pnorm(what + (1/what)*log(vhat/what))
}

pp_doublesaddle_SCC_BN_speedup = function(u,x,muhat){
    nx = dim(x)[2]
    g = x[,1]
    
    x2 = x[,2:nx]
    varx2 = t(muhat*(1-muhat)*x2)%*%x2
    
    g12 = which(g > 0)
    xm = x[g12,]
    
    x2m = xm[,2:nx]
    muhatm = muhat[g12]
    
    varx2m = t(muhatm*(1-muhatm)*x2m)%*%x2m
    var2diff = matrix(0,ncol = nx,nrow = nx)
    var2diff[2:nx,2:nx] = varx2-varx2m
    
    
    x2m = xm[,2:nx]
    var2diff_small = var2diff[2:nx,2:nx]
    
    uoff = u - 0.5
    
    that = optim(rep(0,nx),function(t)k_speedup(t,xm,muhatm,var2diff)-uoff*t[1])$par
    
    what = sign(that[1])*sqrt(2*(-k_speedup(that,xm,muhatm,var2diff)+that[1]*uoff)) 
    
    vhat = 2*sinh(that[1]/2)*sqrt(det(D2k_speedup(that,xm,muhatm,var2diff))/det(D2k_speedup(rep(0,(nx-1)),x2m,muhatm, var2diff_small)))
    c(1 - pnorm(what + (1/what)*log(vhat/what)))
    
}


