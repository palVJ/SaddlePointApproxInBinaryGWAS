# Cumulant generating functions

# For double saddlepoint:
# x is a matrix of covariates under H1 (including intercept)
# and the first column is the genotype

# For efficient score / single saddlepoint:
# x is the adjusted genotype vector

# mu is the vector of estimated means under H0

# Function k: CGF, Returns a scalar
# Function D2k: Hessian of CGF, returns a dxd matrix where d is the number of columns in x

k = function(t,x,mu){
    x = as.matrix(x)
    xmu = mu%*%x
    etx = exp(x%*%t)
    sum(log(1-mu*(1-etx)))-sum(t*xmu)
}


D2k = function(t,x,mu){
    x = as.matrix(x)
    emtx = exp(-x%*%t)
    t(as.numeric((mu*(1-mu)*emtx/(((1-mu)*emtx+mu)^2)))*x)%*%x
}


# Speed-up for double saddlepoint
# xm covariates for individuals with genotype > 0
# mum is the vector of estimated means under H0 for individuals with genotype > 0
# var2diff = var_all - var_m
# is the covariance of the nuisance score for all individuals (pre-computed),
# minus the variance of the nuisance score for individuals with genotype > 0

k_speedup = function(t,xm, mum, var2diff){ 
    xm = as.matrix(xm) 
    xmu = mum%*%xm
    etx = exp(xm%*%t) 
    sum(log(1-mum*(1-etx)))-sum(t*xmu) + 0.5*t(t)%*%var2diff%*%t
}

D2k_speedup = function(t,xm, mum, var2diff){ 
    xm = as.matrix(xm) 
    lent = length(t) 
    emtx = exp(-xm%*%t)
    t(as.numeric((mum*(1-mum)*emtx/(((1-mum)*emtx+mum)^2)))*xm)%*%xm + var2diff 
}