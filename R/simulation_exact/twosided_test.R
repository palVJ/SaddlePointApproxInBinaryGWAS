
# Based on the saddlepoint of the efficient score,
# find two-sided p-value

# u is the observed efficient score
# gadj is the adjusted genotype 
# muhat is a vector of estimated means under H0
# ran is the range of U


singlesaddle2_BN = function(u,gadj,muhat,ran){
  
  if(u > abs(ran[1])){
    # Positive observation and right-tailed test
    
    if(u == ran[2]){
      u = u - 0.5
    }
    
    pu = pp_singlesaddle_BN(u,gadj,muhat)
    return(pu)
  }else if(u > 0){
    # Positive observation and two-tailed test
    uneg = u-sign(u)*ceiling(2*abs(u))
    
    if(uneg == ran[1]){
      uneg = uneg + 0.5
    }
    
    if(u == ran[2]){
      u = u - 0.5
    }
    
    
    pu = pp_singlesaddle_BN(u,gadj,muhat)
    pl = 1 - pp_singlesaddle_BN(uneg,gadj,muhat)
    
    return(pu + pl)
  }else if(u < -1*ran[2]){
    # Negative observation and left-tailed test
    if(u == ran[1]){
      u = u + 0.5
    }
    
    pl = 1 - pp_singlesaddle_BN(u,gadj,muhat)
    return(pl)
  }else{
    # Negative observation and two-tailed test
    upos = u-sign(u)*ceiling(2*abs(u))
    
    if(upos == ran[2]){
      upos = upos - 0.5
    }
    
    if(u == ran[1]){
      u = u + 0.5
    }
    
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


singlesaddle2_SCC_BN_midp = function(u,gadj,muhat,ran){
  
  if((abs(u) <= abs(abs(u)-1)) &  (abs(u)-1 < 0)){
    # observation is as close to 0 as possible, p = 1
    return(1)
  }else if(u > abs(ran[1])){
    # Positive observation and right-tailed test
    pu = pp_singlesaddle_SCC_BN(u,gadj,muhat)
    
    u2 = u + 1
    if(u2 >= ran[2]){
      pu2 = pu
    }else{
      pu2 = pp_singlesaddle_SCC_BN(u2,gadj,muhat)
    }
    return(0.5*(pu+pu2))
  }else if(u > 0){
    # Positive observation and two-tailed test
    uneg = u-sign(u)*ceiling(2*abs(u))
    pu = pp_singlesaddle_SCC_BN(u,gadj,muhat)
    pl = 1 - pp_singlesaddle_SCC_BN(uneg+1,gadj,muhat) # + 1 for survival fn
    
    u2 = u + 1
    if(u2 >= ran[2]){
      pu2 = pu
    }else{
      pu2 = pp_singlesaddle_SCC_BN(u2,gadj,muhat)
    }
    
    uneg2 = uneg - 1
    if(uneg2 <= ran[1]){
      pl2 = pl
    }else{
      pl2 = 1 - pp_singlesaddle_SCC_BN(uneg2+1,gadj,muhat)
    }
    
    return(0.5*(pu+pl+pu2+pl2))
  }else if(u < -1*ran[2]){
    # Negative observation and left-tailed test
    pl = 1 - pp_singlesaddle_SCC_BN(u+1,gadj,muhat) # + 1 for survival fn
    
    u2 = u - 1
    if(u2 <= ran[1]){
      pl2 = pl
    }else{
      pl2 = 1 - pp_singlesaddle_SCC_BN(u2+1,gadj,muhat)
    }
    
    return(0.5*(pl+pl2))
  }else{
    # Negative observation and two-tailed test
    upos = u-sign(u)*ceiling(2*abs(u))
    pu = pp_singlesaddle_SCC_BN(upos,gadj,muhat)
    pl = 1 - pp_singlesaddle_SCC_BN(u+1,gadj,muhat) # + 1 for survival fn
    
    u2 = u - 1
    if(u2 <= ran[1]){
      pl2 = pl
    }else{
      pl2 = 1 - pp_singlesaddle_SCC_BN(u2+1,gadj,muhat)
    }
    
    upos2 = upos + 1
  
    if(upos2 >= ran[2]){
      pu2 = pu
    }else{
      pu2 = pp_singlesaddle_SCC_BN(upos2,gadj,muhat)
    }
    
    
    return(0.5*(pu+pl+pu2+pl2))
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


doublesaddle2_SCC_BN_midp = function(u,x,muhat,ran){
  
  if((abs(u) <= abs(abs(u)-1)) &  (abs(u)-1 < 0)){
    # observation is as close to 0 as possible, p = 1
    return(1)
  }else if(u > abs(ran[1])){
    # Positive observation and right-tailed test
    pu = pp_doublesaddle_SCC_BN(u,x,muhat)
    
    u2 = u + 1
    if(u2 >= ran[2]){
      pu2 = pu
    }else{
      pu2 = pp_doublesaddle_SCC_BN(u2,x,muhat)
    }
    return(0.5*(pu+pu2))
  }else if(u > 0){
    # Positive observation and two-tailed test
    uneg = u-sign(u)*ceiling(2*abs(u))
    pu = pp_doublesaddle_SCC_BN(u,x,muhat)
    pl = 1 - pp_doublesaddle_SCC_BN(uneg+1,x,muhat) # + 1 for survival fn
    
    u2 = u + 1
    if(u2 >= ran[2]){
      pu2 = pu
    }else{
      pu2 = pp_doublesaddle_SCC_BN(u2,x,muhat)
    }
    
    uneg2 = uneg - 1
    if(uneg2 <= ran[1]){
      pl2 = pl
    }else{
      pl2 = 1 - pp_doublesaddle_SCC_BN(uneg2+1,x,muhat)
    }
    
    return(0.5*(pu+pl+pu2+pl2))
  }else if(u < -1*ran[2]){
    # Negative observation and left-tailed test
    pl = 1 - pp_doublesaddle_SCC_BN(u+1,x,muhat) # + 1 for survival fn
    
    u2 = u - 1
    if(u2 <= ran[1]){
      pl2 = pl
    }else{
      pl2 = 1 - pp_doublesaddle_SCC_BN(u2+1,x,muhat)
    }
    
    return(0.5*(pl+pl2))
  }else{
    # Negative observation and two-tailed test
    upos = u-sign(u)*ceiling(2*abs(u))
    pu = pp_doublesaddle_SCC_BN(upos,x,muhat)
    pl = 1 - pp_doublesaddle_SCC_BN(u+1,x,muhat) # + 1 for survival fn
    
    u2 = u - 1
    if(u2 <= ran[1]){
      pl2 = pl
    }else{
      pl2 = 1 - pp_doublesaddle_SCC_BN(u2+1,x,muhat)
    }
    
    upos2 = upos + 1
    
    if(upos2 >= ran[2]){
      pu2 = pu
    }else{
      pu2 = pp_doublesaddle_SCC_BN(upos2,x,muhat)
    }
    
    
    return(0.5*(pu+pl+pu2+pl2))
  }
} 

doublesaddle2_BN = function(u,x,muhat,ran){
  
  if((abs(u) <= abs(abs(u)-1)) &  (abs(u)-1 < 0)){
    # observation is as close to 0 as possible, p = 1
    return(1)
  }else if(u > abs(ran[1])){
    # Positive observation and right-tailed test
    if(u == ran[2]){
      u = u - 0.5
    }
    pu = pp_doublesaddle_BN(u,x,muhat)
    return(pu)
  }else if(u > 0){
    # Positive observation and two-tailed test
    uneg = u-sign(u)*ceiling(2*abs(u))
    
    if(uneg == ran[1]){
      uneg = uneg + 0.5
    }
    
    if(u == ran[2]){
      u = u - 0.5
    }
    
    pu = pp_doublesaddle_BN(u,x,muhat)
    pl = 1 - pp_doublesaddle_BN(uneg,x,muhat) 
    return(pu + pl)
  }else if(u < -1*ran[2]){
    # Negative observation and left-tailed test
    
    if(u == ran[1]){
      u = u + 0.5
    }
    
    pl = 1 - pp_doublesaddle_BN(u,x,muhat) 
    return(pl)
  }else{
    # Negative observation and two-tailed test
    upos = u-sign(u)*ceiling(2*abs(u))
    
    if(upos == ran[2]){
      upos = upos - 0.5
    }
    
    if(u == ran[1]){
      u = u + 0.5
    }
    
    pu = pp_doublesaddle_BN(upos,x,muhat)
    pl = 1 - pp_doublesaddle_BN(u,x,muhat) 
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
    