
# Based on the saddlepoint of the efficient score,
# estimate P(U >= u)

# u is the observed efficient score
# gadj is the adjusted genotype 
# muhat is a vector of estimated means under H0

pp_singlesaddle_LR = function(u,gadj,muhat){
    
    that = optimize(function(t)k(t,gadj,muhat)-u*t,c(-10,10),tol=1e-10)$minimum
    what = sign(that)*sqrt(2*(-k(that,gadj,muhat)+that*u))
    vhat = that*sqrt(D2k(that,gadj,muhat))
    
    c(1 - pnorm(what) - dnorm(what)*(1/what-1/vhat))
}

pp_singlesaddle_BN = function(u,gadj,muhat){
    
    that = optimize(function(t)k(t,gadj,muhat)-u*t,c(-10,10),tol=1e-10)$minimum
    what = sign(that)*sqrt(2*(-k(that,gadj,muhat)+that*u))
    vhat = that*sqrt(D2k(that,gadj,muhat))
    
    c(1 - pnorm(what + (1/what)*log(vhat/what)))
}


pp_singlesaddle_SCC_LR = function(u,gadj,muhat){
    
    uoff = u - 0.5
    
    stilde = optimize(function(t)k(t,gadj,muhat)-uoff*t,c(-10,10),tol=1e-10)$minimum
    vtilde = 2*sinh(stilde/2)*sqrt(D2k(stilde,gadj,muhat))
    wtilde = sign(stilde)*sqrt(2*(-k(stilde,gadj,muhat)+stilde*uoff))
    
    c(1 - pnorm(wtilde) - dnorm(wtilde)*(1/wtilde-1/vtilde))
}

pp_singlesaddle_SCC_BN = function(u,gadj,muhat){
    
    uoff = u - 0.5
    
    stilde = optimize(function(t)k(t,gadj,muhat)-uoff*t,c(-10,10),tol=1e-10)$minimum
    vtilde = 2*sinh(stilde/2)*sqrt(D2k(stilde,gadj,muhat))
    wtilde = sign(stilde)*sqrt(2*(-k(stilde,gadj,muhat)+stilde*uoff))
    
    c(1 - pnorm(wtilde + (1/wtilde)*log(vtilde/wtilde)))
}



