
search_right_ESPA_SCC = function(startRight,uvals,gadj,muhat,ran,pm,alpha){
  
  prob_reject_right = 0
  uid = startRight
  
  tmp = singlesaddle2_SCC_BN(uvals[uid],gadj,muhat,ran)
  while(tmp <= alpha){
    prob_reject_right = sum(pm[uid:length(pm)])
    uid = uid-1
    tmp = singlesaddle2_SCC_BN(uvals[uid],gadj,muhat,ran)
  }
  
  return(prob_reject_right)
}


search_left_ESPA_SCC = function(startLeft,uvals,gadj,muhat,ran,pm,alpha){
  
  prob_reject_left = 0
  uid = startLeft
  
  tmp = singlesaddle2_SCC_BN(uvals[uid],gadj,muhat,ran)
  while(tmp <= alpha){
    prob_reject_left = sum(pm[1:uid])
    uid = uid+1
    tmp = singlesaddle2_SCC_BN(uvals[uid],gadj,muhat,ran)
  }
  
  return(prob_reject_left)
}



############################################################################


search_right_ESPA = function(startRight,uvals,gadj,muhat,ran,pm,alpha){
  
  prob_reject_right = 0
  uid = startRight
  
  tmp = singlesaddle2_BN(uvals[uid],gadj,muhat,ran)
  while(tmp <= alpha){
    prob_reject_right = sum(pm[uid:length(pm)])
    uid = uid-1
    tmp = singlesaddle2_BN(uvals[uid],gadj,muhat,ran)
  }
  
  return(prob_reject_right)
}


search_left_ESPA = function(startLeft,uvals,gadj,muhat,ran,pm,alpha){
  
  prob_reject_left = 0
  uid = startLeft
  
  tmp = singlesaddle2_BN(uvals[uid],gadj,muhat,ran)
  while(tmp <= alpha){
    prob_reject_left = sum(pm[1:uid])
    uid = uid+1
    tmp = singlesaddle2_BN(uvals[uid],gadj,muhat,ran)
  }
  
  return(prob_reject_left)
}


############################################################################



search_right_ESPA_SCC_midp = function(startRight,uvals,gadj,muhat,ran,pm,alpha){
  
  prob_reject_right = 0
  uid = startRight
  
  tmp = singlesaddle2_SCC_BN_midp(uvals[uid],gadj,muhat,ran)
  while(tmp <= alpha){
    prob_reject_right = sum(pm[uid:length(pm)])
    uid = uid-1
    tmp = singlesaddle2_SCC_BN_midp(uvals[uid],gadj,muhat,ran)
  }
  
  return(prob_reject_right)
}


search_left_ESPA_SCC_midp = function(startLeft,uvals,gadj,muhat,ran,pm,alpha){
  
  prob_reject_left = 0
  uid = startLeft
  
  tmp = singlesaddle2_SCC_BN_midp(uvals[uid],gadj,muhat,ran)
  while(tmp <= alpha){
    prob_reject_left = sum(pm[1:uid])
    uid = uid+1
    tmp = singlesaddle2_SCC_BN_midp(uvals[uid],gadj,muhat,ran)
  }
  
  return(prob_reject_left)
}
