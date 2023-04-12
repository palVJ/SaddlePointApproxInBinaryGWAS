

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

