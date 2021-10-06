
# Based on the saddlepoint of the efficient score,
# find two-sided p-value

# u is the observed efficient score
# gadj is the adjusted genotype 
# muhat is a vector of estimated means under H0
# ran is the range of U

# NOTE: cannot handle observation at ran


singlesaddle2_BN = function(u,gadj,muhat,ran){
    
    if(u > abs(ran[1])){
        # Positive observation and right-tailed test
        pu = pp_singlesaddle_BN(u,gadj,muhat)
        return(pu)
    }else if(u > 0){
        # Positive observation and two-tailed test
        uneg = -1*u
        pu = pp_singlesaddle_BN(u,gadj,muhat)
        pl = 1 - pp_singlesaddle_BN(uneg,gadj,muhat)
        return(pu + pl)
    }else if(u < -1*ran[2]){
        # Negative observation and left-tailed test
        pl = 1 - pp_singlesaddle_BN(u,gadj,muhat)
        return(pl)
    }else{
        # Negative observation and two-tailed test
        upos = abs(u)
        pu = pp_singlesaddle_BN(upos,gadj,muhat)
        pl = 1 - pp_singlesaddle_BN(u,gadj,muhat) 
        return(pu + pl)
    }
}


singlesaddle2_SCC_BN = function(u,gadj,muhat,ran){
    
    if((abs(u) <= abs(abs(u)-1)) &  (abs(u)-1 < 0)){
        # observation is as close to 0 as possible, p = 1
        return(1)
    }else if(u > abs(ran[1])){
        # Positive observation and right-tailed test
        pu = pp_singlesaddle_SCC_BN(u,gadj,muhat)
        return(pu)
    }else if(u > 0){
        # Positive observation and two-tailed test
        uneg = u-sign(u)*ceiling(2*abs(u))
        pu = pp_singlesaddle_SCC_BN(u,gadj,muhat)
        pl = 1 - pp_singlesaddle_SCC_BN(uneg+1,gadj,muhat) # + 1 for survival fn
        return(pu + pl)
    }else if(u < -1*ran[2]){
        # Negative observation and left-tailed test
        pl = 1 - pp_singlesaddle_SCC_BN(u+1,gadj,muhat) # + 1 for survival fn
        return(pl)
    }else{
        # Negative observation and two-tailed test
        upos = u-sign(u)*ceiling(2*abs(u))
        pu = pp_singlesaddle_SCC_BN(upos,gadj,muhat)
        pl = 1 - pp_singlesaddle_SCC_BN(u+1,gadj,muhat) # + 1 for survival fn
        return(pu + pl)
    }
} 
    


doublesaddle2_SCC_BN = function(u,x,muhat,ran){
    
    if((abs(u) <= abs(abs(u)-1)) &  (abs(u)-1 < 0)){
        # observation is as close to 0 as possible, p = 1
        return(1)
    }else if(u > abs(ran[1])){
        # Positive observation and right-tailed test
        pu = pp_doublesaddle_SCC_BN(u,x,muhat)
        return(pu)
    }else if(u > 0){
        # Positive observation and two-tailed test
        uneg = u-sign(u)*ceiling(2*abs(u))
        pu = pp_doublesaddle_SCC_BN(u,x,muhat)
        pl = 1 - pp_doublesaddle_SCC_BN(uneg+1,x,muhat) # + 1 for survival fn
        return(pu + pl)
    }else if(u < -1*ran[2]){
        # Negative observation and left-tailed test
        pl = 1 - pp_doublesaddle_SCC_BN(u+1,x,muhat) # + 1 for survival fn
        return(pl)
    }else{
        # Negative observation and two-tailed test
        upos = u-sign(u)*ceiling(2*abs(u))
        pu = pp_doublesaddle_SCC_BN(upos,x,muhat)
        pl = 1 - pp_doublesaddle_SCC_BN(u+1,x,muhat) # + 1 for survival fn
        return(pu + pl)
    }
} 




doublesaddle2_SCC_BN_speedup = function(u,x,muhat,ran){
    
    if((abs(u) <= abs(abs(u)-1)) &  (abs(u)-1 < 0)){
        # observation is as close to 0 as possible, p = 1
        return(1)
    }else if(u > abs(ran[1])){
        # Positive observation and right-tailed test
        pu = pp_doublesaddle_SCC_BN_speedup(u,x,muhat)
        return(pu)
    }else if(u > 0){
        # Positive observation and two-tailed test
        uneg = u-sign(u)*ceiling(2*abs(u))
        pu = pp_doublesaddle_SCC_BN_speedup(upos,x,muhat)
        pl = 1 - pp_doublesaddle_SCC_BN_speedup(uneg+1,x,muhat) # + 1 for survival fn
        return(pu + pl)
    }else if(u < -1*ran[2]){
        # Negative observation and left-tailed test
        pl = 1 - pp_doublesaddle_SCC_BN_speedup(u+1,x,muhat) # + 1 for survival fn
        return(pl)
    }else{
        # Negative observation and two-tailed test
        upos = u-sign(u)*ceiling(2*abs(u))
        pu = pp_doublesaddle_SCC_BN_speedup(upos,x,muhat)
        pl = 1 - pp_doublesaddle_SCC_BN_speedup(u+1,x,muhat) # + 1 for survival fn
        return(pu + pl)
    }
} 
    