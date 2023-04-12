
search_right_DSPA_SCC = function(startRight,uvals,x,muhat,ran,pm,alpha){
  
  prob_reject_right = 0
  uid = startRight
  
  tmp = doublesaddle2_SCC_BN(uvals[uid],x,muhat,ran)
  while(tmp < alpha){
    prob_reject_right = sum(pm[uid:length(pm)])
    uid = uid-1
    tmp = doublesaddle2_SCC_BN(uvals[uid],x,muhat,ran)
  }
  
  return(prob_reject_right)
}


search_left_DSPA_SCC = function(startLeft,uvals,x,muhat,ran,pm,alpha){
  
  prob_reject_left = 0
  uid = startLeft
  
  tmp = doublesaddle2_SCC_BN(uvals[uid],x,muhat,ran)
  while(tmp < alpha){
    prob_reject_left = sum(pm[1:uid])
    uid = uid+1
    tmp = doublesaddle2_SCC_BN(uvals[uid],x,muhat,ran)
  }
  
  return(prob_reject_left)
}



########################################################################################


search_right_DSPA = function(startRight,uvals,x,muhat,ran,pm,alpha){
  
  prob_reject_right = 0
  uid = startRight
  
  tmp = doublesaddle2_BN(uvals[uid],x,muhat,ran)
  while(tmp < alpha){
    prob_reject_right = sum(pm[uid:length(pm)])
    uid = uid-1
    tmp = doublesaddle2_BN(uvals[uid],x,muhat,ran)
  }
  
  return(prob_reject_right)
}


search_left_DSPA = function(startLeft,uvals,x,muhat,ran,pm,alpha){
  
  prob_reject_left = 0
  uid = startLeft
  
  tmp = doublesaddle2_BN(uvals[uid],x,muhat,ran)
  while(tmp < alpha){
    prob_reject_left = sum(pm[1:uid])
    uid = uid+1
    tmp = doublesaddle2_BN(uvals[uid],x,muhat,ran)
  }
  
  return(prob_reject_left)
}



########################################################################################


search_right_DSPA_SCC_midp = function(startRight,uvals,x,muhat,ran,pm,alpha){
  
  prob_reject_right = 0
  uid = startRight
  
  tmp = doublesaddle2_SCC_BN_midp(uvals[uid],x,muhat,ran)
  while(tmp < alpha){
    prob_reject_right = sum(pm[uid:length(pm)])
    uid = uid-1
    tmp = doublesaddle2_SCC_BN_midp(uvals[uid],x,muhat,ran)
  }
  
  return(prob_reject_right)
}


search_left_DSPA_SCC_midp = function(startLeft,uvals,x,muhat,ran,pm,alpha){
  
  prob_reject_left = 0
  uid = startLeft
  
  tmp = doublesaddle2_SCC_BN_midp(uvals[uid],x,muhat,ran)
  while(tmp < alpha){
    prob_reject_left = sum(pm[1:uid])
    uid = uid+1
    tmp = doublesaddle2_SCC_BN_midp(uvals[uid],x,muhat,ran)
  }
  
  return(prob_reject_left)
}