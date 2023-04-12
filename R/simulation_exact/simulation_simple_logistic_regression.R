
# functions
source("CGFs.R")          # cumulant generating functions
source("singlesadle.R")   # single sadlepoint logistic regression with efficient score
source("doublesadle.R")   # double sadlepoint logistic regression 
source("twosided_test.R") # two-sided p-values from single or double sadlepoint
source("exact_test_intercept.R") # exact test 
source("find_significance_doublesadle.R") # find critical regions and 
source("find_significance_singlesadle.R")

###############################################################################

alpha = 5e-8

n = 10000
maf = 0.1

n2 = n*maf^2
n1 = n*2*maf*(1-maf)
n0 = n*(1-maf)^2

g = c(rep(0,n0),rep(1,n1),rep(2,n2))
gtilde = g - mean(g)
x = cbind(g,rep(1,n))




####################################################################################


# exact point probabilities and cumulative probabilities

ncases = 1:(n-1)
u_list = list()
pmf_list = list()
right_tail_list = list()
left_tail_list = list()

for(nc in ncases){
  if(nc %% 100 == 0){
    print(nc)
  }
  mu = nc/n    # muhat for dette antall kasus
  muhat = rep(mu,n) # vektor som brukes i tester
  ran = c(-sum(g)*mu, (1-mu)*sum(g)) # range
  uvals = seq(ran[1],ran[2],1) # alle u-verdier i range
  pm = pmf_condscore(uvals,n0,n1,n2,nc,mu) # eksakte sannsynligheter
  
  pm_pos = which(pm > 0)
  ran_true = c(uvals[pm_pos[1]],uvals[pm_pos[length(pm_pos)]])
  uvals = seq(ran_true[1],ran_true[2],1)
  ran = ran_true
  
  pm = pm[pm_pos[1]:pm_pos[length(pm_pos)]]

  u_list[[nc]] = uvals
  pmf_list[[nc]] = pm
  
  tmp = length(pm)
  right_tail_list[[nc]] = cumsum(pm[tmp:1])[tmp:1]
  left_tail_list[[nc]] = cumsum(pm)
}


####################################################################################


# rejection region and tail probabilities for two-sided exact test

critical_right2 = c()
critical_left2 = c()
prob_reject2 = c()

critical_right2_index = c()
critical_left2_index = c()

for(nc in ncases){
  if(nc %% 100 == 0){
    print(nc)
  }
  rtailprob = right_tail_list[[nc]]
  ltailprob = left_tail_list[[nc]]
  uvals = u_list[[nc]]
  len = length(uvals)
  
  critical_left2[nc] = NA
  critical_right2[nc] = NA
  critical_right2_index[nc] = NA
  critical_left2_index[nc] = NA
  
  
  prob_reject2[nc] = 0
  
  if(ltailprob[1]>alpha){ # da gjøres bare - hvis mulig - test i høyre hale
    tmp = which(rtailprob<alpha)
    if(length(tmp)>0){
      critical_right2[nc] = uvals[min(tmp)]
      prob_reject2[nc] = rtailprob[min(tmp)]
      critical_right2_index[nc] = min(tmp)
    }
  }else if(rtailprob[len]>alpha){ # left tail only
    tmp = which(ltailprob<alpha)
    if(length(tmp)>0){
      critical_left2[nc] = uvals[max(tmp)]
      prob_reject2[nc] = ltailprob[max(tmp)]
      critical_left2_index[nc] = max(tmp)
    }
  }else{ # two-sided
    upos = uvals[len]
    uneg = uvals[1]
    uopp = upos - ceiling(2*upos)
    
    c = len
    
    while( (uopp < uneg) ){
      c = c - 1
      upos = uvals[c]
      uopp = upos - ceiling(2*upos)
    }
    
    r = which(round(uvals,2) == round(uopp,2))
    
    if((ltailprob[r] + rtailprob[c] > alpha) & (uvals[len] > abs(uvals[1])) ){ 
      tmp = which(rtailprob<alpha)
      critical_right2[nc] = uvals[min(tmp[tmp>c])] # her forkastes for høyre-observasjoner, ensidig
      prob_reject_right = rtailprob[min(tmp[tmp>c])]
      critical_right2_index[nc] = min(tmp[tmp>c])
      
    }else if((ltailprob[r] + rtailprob[c] > alpha) & (uvals[len] <= abs(uvals[1])) ){
      prob_reject_right = 0
      critical_right2_index[nc] = NA
    }else{ 
      probtmp = rtailprob[c] + ltailprob[r]
      
      while(probtmp < alpha){
        c = c - 1
        r = r + 1
        probtmp = rtailprob[c] + ltailprob[r]
      }
      c = c + 1
      r = r - 1
      
      critical_right2[nc] = uvals[c]
      ind_right = c
      
      prob_reject_right = rtailprob[c]
      
      critical_right2_index[nc] = c
    }
  
    upos = uvals[len]
    uneg = uvals[1]
    uopp = uneg + ceiling(2*abs(uneg))
    
    r = 1
    
    while( (uopp > upos) ){
      r = r + 1
      uneg = uvals[r]
      uopp = uneg + ceiling(2*abs(uneg))
    }
    
    c = which(round(uvals,2) == round(uopp,2))
    
    if((ltailprob[r] + rtailprob[c] > alpha)& (uvals[len] < abs(uvals[1]))){ 
      
      tmp = which(ltailprob<alpha)
      
      critical_left2[nc] = uvals[max(tmp[tmp<r])] 
      prob_reject_left = ltailprob[max(tmp[tmp<r])]
      critical_left2_index[nc] = max(tmp[tmp<r])
      
    }else if((ltailprob[r] + rtailprob[c] > alpha) & (uvals[len] >= abs(uvals[1])) ){
      prob_reject_left = 0
    }else{ 
      probtmp = rtailprob[c] + ltailprob[r]
      while(probtmp < alpha){
        c = c - 1
        r = r + 1
        probtmp = rtailprob[c] + ltailprob[r]
      }
      c = c + 1 
      r = r - 1
      
      critical_left2[nc] = uvals[r]
      critical_left2_index[nc] = r
      ind_left = r
      
      prob_reject_left = ltailprob[r]
    }
    prob_reject2[nc] = prob_reject_left + prob_reject_right
  }
  
}



# Same using mid-p 

critical_right2_midp = c()
critical_left2_midp = c()
prob_reject2_midp = c()
critical_right2_midp_index = c()
critical_left2_midp_index = c()

for(nc in ncases){
  if(nc %% 100 == 0){
    print(nc)
  }
  rtailprob_exact = right_tail_list[[nc]]
  
  rtailprob = rtailprob_exact
  rtailprob[1:(length(rtailprob_exact)-1)] = 0.5*rtailprob_exact[1:(length(rtailprob_exact)-1)] + 0.5*rtailprob_exact[2:length(rtailprob_exact)]
  
  
  ltailprob_exact = left_tail_list[[nc]]
  ltailprob = ltailprob_exact
  
  ltailprob[2:length(ltailprob_exact)] = 0.5*ltailprob_exact[2:length(ltailprob_exact)] + 0.5*ltailprob_exact[1:(length(ltailprob_exact)-1)]

  
  uvals = u_list[[nc]]
  len = length(uvals)
  
  
  prob_reject2_midp[nc] = 0
  
  if(ltailprob[1]>alpha){ 
    tmp = which(rtailprob<alpha)
    if(length(tmp)>0){
      critical_right2_midp[nc] = uvals[min(tmp)]
      prob_reject2_midp[nc] = rtailprob_exact[min(tmp)]
      critical_right2_midp_index[nc] = min(tmp)
    }
  }else if(rtailprob[len]>alpha){ 
    tmp = which(ltailprob<alpha)
    if(length(tmp)>0){
      critical_left2_midp[nc] = uvals[max(tmp)]
      prob_reject2_midp[nc] = ltailprob_exact[max(tmp)]
      critical_left2_midp_index[nc] = max(tmp)
    }
  }else{ 
    upos = uvals[len]
    uneg = uvals[1]
    uopp = upos - ceiling(2*upos)
    
    c = len
    
    while( (uopp < uneg) ){
      c = c - 1
      upos = uvals[c]
      uopp = upos - ceiling(2*upos)
    }
    
    r = which(round(uvals,2) == round(uopp,2))
    
    if((ltailprob[r] + rtailprob[c] > alpha) & (uvals[len] > abs(uvals[1])) ){ 
      
      
      tmp = which(rtailprob<alpha)
      critical_right2_midp[nc] = uvals[min(tmp[tmp>c])] 
      prob_reject_right = rtailprob_exact[min(tmp[tmp>c])]
      critical_right2_midp_index[nc] = min(tmp[tmp>c])
      
    }else if((ltailprob[r] + rtailprob[c] > alpha) & (uvals[len] <= abs(uvals[1])) ){
      prob_reject_right = 0
      critical_right2_midp_index[nc] = NA
    }else{ 
      probtmp = rtailprob[c] + ltailprob[r]
      
      while(probtmp < alpha){
        c = c - 1
        r = r + 1
        probtmp = rtailprob[c] + ltailprob[r]
      }
      c = c + 1 
      r = r - 1
      
      critical_right2_midp[nc] = uvals[c]
      ind_right = c
      
      prob_reject_right = rtailprob_exact[c]
      
      critical_right2_midp_index[nc] = c
    }

    upos = uvals[len]
    uneg = uvals[1]
    uopp = uneg + ceiling(2*abs(uneg))
    
    r = 1
    
    while( (uopp > upos) ){
      r = r + 1
      uneg = uvals[r]
      uopp = uneg + ceiling(2*abs(uneg))
    }
    
    c = which(round(uvals,2) == round(uopp,2))
    
    if((ltailprob[r] + rtailprob[c] > alpha)& (uvals[len] < abs(uvals[1]))){ 
      
      
      tmp = which(ltailprob<alpha)
      
      critical_left2_midp[nc] = uvals[max(tmp[tmp<r])] 
      prob_reject_left = ltailprob_exact[max(tmp[tmp<r])]
      critical_left2_midp_index[nc] = max(tmp[tmp<r])
      
    }else if((ltailprob[r] + rtailprob[c] > alpha) & (uvals[len] >= abs(uvals[1])) ){
      prob_reject_left = 0
    }else{ 
      probtmp = rtailprob[c] + ltailprob[r]
      while(probtmp < alpha){
        c = c - 1
        r = r + 1
        probtmp = rtailprob[c] + ltailprob[r]
      }
      c = c + 1 
      r = r - 1
      
      critical_left2_midp[nc] = uvals[r]
      critical_left2_midp_index[nc] = r
      ind_left = r
      
      prob_reject_left = ltailprob_exact[r]
    }
    prob_reject2_midp[nc] = prob_reject_left + prob_reject_right
  }
  
}

####################################################################################

# Single sadlepoint, efficient score

critical_right2_ESPA = c()
critical_left2_ESPA = c()
prob_reject2_ESPA = c()

for(nc in ncases){
    if(nc %% 100 == 0){
        print(nc)
    }
  pm = pmf_list[[nc]]
  uvals = u_list[[nc]]
  len = length(uvals)
  
  prob_reject_right = 0
  prob_reject_left = 0
  
  mu = nc/n   
  muhat = rep(mu,n) 
  ran = c(min(uvals), max(uvals)) 
  
  critical_left2_ESPA[nc] = NA
  critical_right2_ESPA[nc] = NA
  
  exactRight = critical_right2_index[nc]
  exactLeft = critical_left2_index[nc]
  
  if(is.na(exactRight) & is.na(exactLeft)){

    prob_reject_right = search_right_ESPA_SCC(len,uvals,gtilde,muhat,ran,pm,alpha)
    
    prob_reject_left = search_left_ESPA_SCC(1,uvals,gtilde,muhat,ran,pm,alpha)
    
  }else if(is.na(exactLeft)){
    startid = min(exactRight+2,len)
    prob_reject_right = search_right_ESPA_SCC(startid,uvals,gtilde,muhat,ran,pm,alpha)
    
    upos = uvals[exactRight]; uneg = uvals[1]; uopp = upos - ceiling(2*upos)
    if(uneg-2 <= uopp){
      prob_reject_left = search_left_ESPA_SCC(1,uvals,gtilde,muhat,ran,pm,alpha)
    }
  }else if (is.na(exactRight)){
    startid = max(exactLeft-3,1)
    prob_reject_left = search_left_ESPA_SCC(startid,uvals,gtilde,muhat,ran,pm,alpha)
    
    upos = uvals[len]; uneg = uvals[exactLeft]; uopp = uneg + ceiling(2*abs(uneg))
    
    if(upos+3 >= uopp){
      prob_reject_right = search_right_ESPA_SCC(len,uvals,gtilde,muhat,ran,pm,alpha)
    }
    
  }else{
    startid = min(exactRight+2,len)
    prob_reject_right = search_right_ESPA_SCC(startid,uvals,gtilde,muhat,ran,pm,alpha)
    
    startid = max(exactLeft-2,1)
    prob_reject_left = search_left_ESPA_SCC(startid,uvals,gtilde,muhat,ran,pm,alpha)
  }
  
  prob_reject2_ESPA[nc] = prob_reject_left + prob_reject_right
}

####################################################################################

prob_reject2_ESPAnaive = c()

for(nc in ncases){
    if(nc %% 100 == 0){
        print(nc)
    }
  pm = pmf_list[[nc]]
  uvals = u_list[[nc]]
  len = length(uvals)
  
  prob_reject_right = 0
  prob_reject_left = 0
  
  mu = nc/n    
  muhat = rep(mu,n) 
  ran = c(min(uvals), max(uvals)) 
  
  exactRight = critical_right2_index[nc]
  exactLeft = critical_left2_index[nc]
  
  if(is.na(exactRight) & is.na(exactLeft)){

    prob_reject_right = search_right_ESPA(len,uvals,gtilde,muhat,ran,pm,alpha)
    
    prob_reject_left = search_left_ESPA(1,uvals,gtilde,muhat,ran,pm,alpha)
    
  }else if(is.na(exactLeft)){
    startid = min(exactRight+2,len)
    prob_reject_right = search_right_ESPA(startid,uvals,gtilde,muhat,ran,pm,alpha)
    
    upos = uvals[exactRight]; uneg = uvals[1]; uopp = upos - ceiling(2*upos)
    if(uneg-2 <= uopp){

      prob_reject_left = search_left_ESPA(1,uvals,gtilde,muhat,ran,pm,alpha)
    }
  }else if (is.na(exactRight)){
    startid = max(exactLeft-2,1)
    prob_reject_left = search_left_ESPA(startid,uvals,gtilde,muhat,ran,pm,alpha)
    
    upos = uvals[len]; uneg = uvals[exactLeft]; uopp = uneg + ceiling(2*abs(uneg))
    
    if(upos+2 >= uopp){

      prob_reject_right = search_right_ESPA(len,uvals,gtilde,muhat,ran,pm,alpha)
    }
    
  }else{
    startid = min(exactRight+2,len)
    prob_reject_right = search_right_ESPA(startid,uvals,gtilde,muhat,ran,pm,alpha)
    
    startid = max(exactLeft-2,1)
    prob_reject_left = search_left_ESPA(startid,uvals,gtilde,muhat,ran,pm,alpha)
  }
  
  prob_reject2_ESPAnaive[nc] = prob_reject_left + prob_reject_right
}

####################################################################################
# Single sadelpoint, efficients score, cont. correction and mid-p
prob_reject2_ESPA_midp = c()

for(nc in ncases){
    if(nc %% 100 == 0){
        print(nc)
    }
  pm = pmf_list[[nc]]
  uvals = u_list[[nc]]
  len = length(uvals)
  
  prob_reject_right = 0
  prob_reject_left = 0
  
  mu = nc/n    
  muhat = rep(mu,n) 
  ran = c(min(uvals), max(uvals)) 
  
  exactRight = critical_right2_index[nc]
  exactLeft = critical_left2_index[nc]
  
  if(is.na(exactRight) & is.na(exactLeft)){

    prob_reject_right = search_right_ESPA_SCC_midp(len,uvals,gtilde,muhat,ran,pm,alpha)
    
    prob_reject_left = search_left_ESPA_SCC_midp(1,uvals,gtilde,muhat,ran,pm,alpha)
    
  }else if(is.na(exactLeft)){
    startid = min(exactRight+2,len)
    prob_reject_right = search_right_ESPA_SCC_midp(startid,uvals,gtilde,muhat,ran,pm,alpha)
    
    upos = uvals[exactRight]; uneg = uvals[1]; uopp = upos - ceiling(2*upos)
    if(uneg-2 <= uopp){

      prob_reject_left = search_left_ESPA_SCC_midp(1,uvals,gtilde,muhat,ran,pm,alpha)
    }
  }else if (is.na(exactRight)){
    startid = max(exactLeft-2,1)
    prob_reject_left = search_left_ESPA_SCC_midp(startid,uvals,gtilde,muhat,ran,pm,alpha)
    
    upos = uvals[len]; uneg = uvals[exactLeft]; uopp = uneg + ceiling(2*abs(uneg))
    
    if(upos+2 >= uopp){

      prob_reject_right = search_right_ESPA_SCC_midp(len,uvals,gtilde,muhat,ran,pm,alpha)
    }
    
  }else{
    startid = min(exactRight+2,len)
    prob_reject_right = search_right_ESPA_SCC_midp(startid,uvals,gtilde,muhat,ran,pm,alpha)
    
    startid = max(exactLeft-2,1)
    prob_reject_left = search_left_ESPA_SCC_midp(startid,uvals,gtilde,muhat,ran,pm,alpha)
  }
  
  prob_reject2_ESPA_midp[nc] = prob_reject_left + prob_reject_right
}


####################################################################################

# Double sadlepoint with cont. correction

critical_right2_DSPA = c()
critical_left2_DSPA = c()
prob_reject2_DSPA = c()

for(nc in ncases){
  print(nc)
  pm = pmf_list[[nc]]
  uvals = u_list[[nc]]
  len = length(uvals)

  prob_reject_right = 0
  prob_reject_left = 0

  mu = nc/n    
  muhat = rep(mu,n) 
  ran = c(min(uvals), max(uvals)) 

  critical_left2_DSPA[nc] = NA
  critical_right2_DSPA[nc] = NA

  exactRight = critical_right2_index[nc]
  exactLeft = critical_left2_index[nc]

  if(is.na(exactRight) & is.na(exactLeft)){
    prob_reject_right = search_right_DSPA_SCC(len,uvals,x,muhat,ran,pm,alpha)
    prob_reject_left = search_left_DSPA_SCC(1,uvals,x,muhat,ran,pm,alpha)

  }else if(is.na(exactLeft)){
    startid = min(exactRight+2,len)
    prob_reject_right = search_right_DSPA_SCC(startid,uvals,x,muhat,ran,pm,alpha)

    upos = uvals[exactRight]; uneg = uvals[1]; uopp = upos - ceiling(2*upos)
    if(uneg-2 <= uopp){
      prob_reject_left = search_left_DSPA_SCC(1,uvals,x,muhat,ran,pm,alpha)
    }
  }else if (is.na(exactRight)){
    startid = max(exactLeft-2,1)
    prob_reject_left = search_left_DSPA_SCC(startid,uvals,x,muhat,ran,pm,alpha)

    upos = uvals[len]; uneg = uvals[exactLeft]; uopp = uneg + ceiling(2*abs(uneg))

    if(upos+2 >= uopp){
      prob_reject_right = search_right_DSPA_SCC(len,uvals,x,muhat,ran,pm,alpha)
    }

  }else{
    startid = min(exactRight+2,len)
    prob_reject_right = search_right_DSPA_SCC(startid,uvals,x,muhat,ran,pm,alpha)

    startid = max(exactLeft-2,1)
    prob_reject_left = search_left_DSPA_SCC(startid,uvals,x,muhat,ran,pm,alpha)
  }

  prob_reject2_DSPA[nc] = prob_reject_left + prob_reject_right
}


####################################################################################

# Double sadle with cont. correction and mid-p


prob_reject2_DSPA_midp = c()

for(nc in ncases){
  print(nc)
  pm = pmf_list[[nc]]
  uvals = u_list[[nc]]
  len = length(uvals)
  
  prob_reject_right = 0
  prob_reject_left = 0
  
  mu = nc/n    
  muhat = rep(mu,n) 
  ran = c(min(uvals), max(uvals)) 
  
  exactRight = critical_right2_index[nc]
  exactLeft = critical_left2_index[nc]
  
  if(is.na(exactRight) & is.na(exactLeft)){

    prob_reject_right = search_right_DSPA_SCC_midp(len,uvals,x,muhat,ran,pm,alpha)
    
    prob_reject_left = search_left_DSPA_SCC_midp(1,uvals,x,muhat,ran,pm,alpha)
    
  }else if(is.na(exactLeft)){
    startid = min(exactRight+2,len)
    prob_reject_right = search_right_DSPA_SCC_midp(startid,uvals,x,muhat,ran,pm,alpha)
    
    upos = uvals[exactRight]; uneg = uvals[1]; uopp = upos - ceiling(2*upos)
    if(uneg-2 <= uopp){

      prob_reject_left = search_left_DSPA_SCC_midp(1,uvals,x,muhat,ran,pm,alpha)
    }
  }else if (is.na(exactRight)){
    startid = max(exactLeft-2,1)
    prob_reject_left = search_left_DSPA_SCC_midp(startid,uvals,x,muhat,ran,pm,alpha)
    
    upos = uvals[len]; uneg = uvals[exactLeft]; uopp = uneg + ceiling(2*abs(uneg))
    
    if(upos+2 >= uopp){

      prob_reject_right = search_right_DSPA_SCC_midp(len,uvals,x,muhat,ran,pm,alpha)
    }
    
  }else{
    startid = min(exactRight+2,len)
    prob_reject_right = search_right_DSPA_SCC_midp(startid,uvals,x,muhat,ran,pm,alpha)
    
    startid = max(exactLeft-2,1)
    prob_reject_left = search_left_DSPA_SCC_midp(startid,uvals,x,muhat,ran,pm,alpha)
  }
  
  prob_reject2_DSPA_midp[nc] = prob_reject_left + prob_reject_right
}

####################################################################################


# Normalapprox

prob_reject2_norm = c()
for(nc in ncases){
  print(nc)
  pm = pmf_list[[nc]]
  uvals = u_list[[nc]]
  len = length(uvals)
  
  ran = c(min(uvals),max(uvals))
  
  mu = nc/n
  
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
    prob_reject2_norm[nc] = sum(pm[1:indneg]) + sum(pm[indpos:length(pm)])
  }else if (crit_norm_inv<ran[1] & crit_norm <= ran[2]){
    indpos = which(round(uvals,2) == round(crit_norm,2))
    prob_reject2_norm[nc] = sum(pm[indpos:length(pm)])
  }else if (crit_norm_inv>=ran[1] & crit_norm > ran[2]){
    indneg = which(round(uvals,2) == round(crit_norm_inv,2))
    prob_reject2_norm[nc] = sum(pm[1:indneg])
  }else{
    prob_reject2_norm[nc] = 0
  }
}



# PLOT

library(ggplot2)
library(reshape2)
library(ggpubr)


yn1 = 1:(n-1)  

mumin = 0.01
mumax = 0.5 
mus = seq(mumin,mumax,0.001)

t1error = c()
probab_invalid = c()

t1error_DSPACC = c()
probab_invalid_DSPACC = c()

t1error_DSPACC_midp = c()
probab_invalid_DSPACC_midp = c()

t1error_DSPA = c()
probab_invalid_DSPA = c()


t1error_ESPACC = c()
probab_invalid_ESPACC = c()

t1error_ESPACC_midp = c()
probab_invalid_ESPACC_midp = c()

t1error_ESPA = c()
probab_invalid_ESPA = c()

t1error_norm = c()
probab_invalid_norm = c()


muno = 0

for(mu in mus){
  muno = muno + 1
  pn1 = dbinom(yn1,n,mu) 
  
  t1error[muno] = sum(exact*pn1)
  probab_invalid[muno] = sum((exact>alpha)*pn1)
  
  t1error_DSPACC[muno] = sum(DSPACC*pn1)
  probab_invalid_DSPACC[muno] = sum((DSPACC>alpha)*pn1)
  
  t1error_DSPACC_midp[muno] = sum(DSPACCmidp*pn1)
  probab_invalid_DSPACC_midp[muno] = sum((DSPACCmidp>alpha)*pn1)
  
  t1error_DSPA[muno] = sum(DSPA*pn1)
  probab_invalid_DSPA[muno] = sum((DSPA>alpha)*pn1)
  
  t1error_ESPACC[muno] = sum(ESPACC*pn1)
  probab_invalid_ESPACC[muno] = sum((ESPACC>alpha)*pn1)
  
  t1error_ESPACC_midp[muno] = sum(ESPACCmidp*pn1)
  probab_invalid_ESPACC_midp[muno] = sum((ESPACCmidp>alpha)*pn1)
  
  t1error_ESPA[muno] = sum(ESPA*pn1)
  probab_invalid_ESPA[muno] = sum((ESPA>alpha)*pn1)
  
  
  t1error_norm[muno] = sum(norm*pn1)
  probab_invalid_norm[muno] = sum((norm>alpha)*pn1)
  
}



cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


plot(x = mus, y = t1error,type = "l", ylim = c(alpha*0.5,alpha+alpha*0.5), lwd = 4, col = "black",yaxt = "n",
     xaxt = "n", xlab = expression(mu), ylab = "P(type I error)")
axis(1, at = seq(0,1,0.1))
axis(2, at = seq(alpha*0.5,alpha+alpha*0.5,alpha/2),labels = c("2.5e-8", "5.0e-8", "7.5e-8"))
abline(h = alpha,lty = 1, col = "grey20")
lines(x = mus, y = t1error,lty=1, lwd = 4, col = "black")
lines(x = mus, y = t1error_DSPACC, col = cbp1[2], lty = 1, lwd = 2)
lines(x = mus, y = t1error_ESPACC,lty=2, col = cbp1[3] ,lwd = 2)
lines(x = mus, y = t1error_midp,lty=1, lwd = 4, col = cbp1[1])
lines(x = mus, y = t1error_DSPA, col = cbp1[7], lty = 1, lwd = 2)
lines(x = mus, y = t1error_ESPA,lty=2, col = cbp1[6] ,lwd = 2)
lines(x = mus, y = t1error_norm,lty=3, col = cbp1[4] ,lwd = 2)

legend("topright",  # Position
       legend = c("Exact", "DSPA-CC", "ESPA-CC", "Normal", "Exact mid-p", "DSPA", "ESPA"),
       lty = c(1,1,2,3,1,1,2),
       col = c("black",cbp1[2],cbp1[3],cbp1[4],cbp1[1],cbp1[7],cbp1[6]),
       lwd = c(4,2,2,2,4,2,2),
       cex = 0.8,
       bty = "n",
       ncol = 2)



plot(x = mus, y = probab_invalid,type = "l", ylim = c(0,1), lwd = 4, col = "black",yaxt = "n",
     xaxt = "n", xlab = expression(mu),ylab = "")
axis(1, at = seq(0,1,0.1))
axis(2, at = seq(0,1,0.1))
lines(x = mus, y = probab_invalid_DSPACC, col = cbp1[2], lty = 1, lwd = 2)
lines(x = mus, y = probab_invalid_ESPACC,lty=2, col = cbp1[3] ,lwd = 2)

lines(x = mus, y = probab_invalid_midp, col = cbp1[1], lty = 1, lwd = 4)
lines(x = mus, y = probab_invalid_DSPA, col = cbp1[7], lty = 1, lwd = 2)
lines(x = mus, y = probab_invalid_ESPA,lty=2, col = cbp1[6] ,lwd = 2)

lines(x = mus, y = probab_invalid_norm,lty=3, col = cbp1[4] ,lwd = 2)

title(ylab="P(conditionally invalid test)", line=2, cex.lab=1.2)
