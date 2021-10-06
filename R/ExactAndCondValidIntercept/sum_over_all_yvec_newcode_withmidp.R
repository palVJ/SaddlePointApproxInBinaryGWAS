

start = Sys.time()

library(ggplot2)
library(reshape2)
library(ggpubr)
library(viridis)

source("~/kode_intercept copy/CGFs.R")
source("~/kode_intercept copy/singlesadle.R")
source("~/kode_intercept copy/doublesadle.R")
source("~/kode_intercept copy/twosided_test.R")
source("~/kode_intercept copy/exact_test_intercept.R")
source("~/kode_intercept copy/exact_t1error.R")



midp = function(uid,uvals,ran,pm){
    u = uvals[uid]
    
    if(u > abs(ran[1])){
        # Positive observation and right-tailed test
        
        if(uid == length(pm)){
            pu = pm[length(pm)]
        }else{
            pu = 0.5*sum(pm[uid:length(pm)]) + 0.5*sum(pm[(uid+1):length(pm)])
        }
        return(pu)
    }else if(u > 0){
        # Positive observation and two-tailed test
        uneg = u-sign(u)*ceiling(2*abs(u))
        unegid = which(round(uvals,2) == round(uneg,2))
        
        if(uid == length(pm)){
            pu = pm[length(pm)]
        }else{
            pu = 0.5*sum(pm[uid:length(pm)]) + 0.5*sum(pm[(uid+1):length(pm)])
        }
        
        if(unegid == 1){
            pl = pm[1]
        }else{
            pl = 0.5*sum(pm[1:unegid]) + 0.5*sum(pm[1:(unegid-1)])
        }
        return(pu + pl)
    }else if(u < -1*ran[2]){
        # Negative observation and left-tailed test
        
        if(uid == 1){
            pl = pm[1]
        }else{
            pl = 0.5*sum(pm[1:uid]) + 0.5*sum(pm[1:(uid-1)])
        }

        return(pl)
    }else{
        # Negative observation and two-tailed test
        upos = u-sign(u)*ceiling(2*abs(u))
        uposid = which(round(uvals,2) == round(upos,2))
        
        if(uposid == length(pm)){
            pu = pm[length(pm)]
        }else{
            pu = 0.5*sum(pm[uposid:length(pm)]) + 0.5*sum(pm[(uposid+1):length(pm)])
        }
        
        if(uid == 1){
            pl = pm[1]
        }else{
            pl = 0.5*sum(pm[1:uid]) + 0.5*sum(pm[1:(uid-1)])
        }
        return(pu + pl)
    }
}



# Regner ut sannsynlighet for type-1 feil i intercept-modell


# Velg signifikansnivå,
# n (sample size)
# velg komposisjon av g-vektor (n0, n1, n2)

# Merk: veldig lave signifikansnivå med liten n 
# kan føre til feilmeldinger i kode, fordi ingen
# av u-verdiene i verdimengen gir liten nok sannsynlighet
n = 1000
n2 = 0; n1 = 20; n0 = 980
g = c(rep(0,n0),rep(1,n1),rep(2,n2))
gtilde = g - mean(g)
x = cbind(g,rep(1,n))

for(r in 1:2){

if(r == 1){
    alpha = 0.05    
}else{
    alpha = 5e-05
}

ncases = 1:(n-1)


# I disse vektorene skal det lagres
# eksakt sannsynlighet for type-1 feil for de ulike
# test-metodene

signExact = c()
signSPA = c()
signSPAcc = c()
signDSPAcc = c()
signNorm = c()
signMidP = c()

c = 0
for(nc in ncases){
    c = c + 1
    
    print(nc)
    
    mu = nc/n    # muhat for dette antall kasus
    
    muhat = rep(mu,n) # vektor som brukes i tester
    
    ran = round(c(-sum(g)*mu, (1-mu)*sum(g)),2) # range
    uvals = round(seq(ran[1],ran[2],1),2) # alle u-verdier i range
    
    pm = pmf_condscore(uvals,n0,n1,n2,nc,mu) # eksakte sannsynligheter
    #Correct the range as ran is lower and upper bounds, not exact bounds:
    #We can find the true range by looking at which uvals that have zero probability according to exact distribution
    gtz = which(pm > 0)
    ran_true = c(uvals[gtz[1]],uvals[gtz[length(gtz)]])
    uvals_true = round(seq(ran_true[1],ran_true[2],1),2)
    
    # En funksjon som returnerer sanns. for type-1 feil til 
    # eksakt test, i tillegg til øvre t1[2] og nedre t1[3]
    # kritiske verdi (NA hvis ikke øvre/nedre fører til forkastning)
    
    t1 = find_t1error(ran,uvals,pm,alpha)  
    
    signExact[c] = t1[1] # eksakt test type-1 feil
    crit = t1[2] # kritisk verdi i høyre hale fra eksakt test
    
    
    ############################################################    
    if( !is.na(crit) & (((uvals[crit]-3) > abs(ran[1])) | (uvals[crit]-3) < 0)){
        # hvis TRUE:
        # Fordelingen til testobservator er right-skew
        # og vi startet søk etter (approksimerte) kritiske 
        # så langt ut i høyre hale at bare høyre-sidig test er vurdert
    ############################################################ 
        
        
        #################
        # SPA
        #################
        uid = max(crit-3,which(uvals>0)[1]) 
        
        # uid er der vi starter vi søket vårt, tre verdier lenger inn mot
        # midten enn eksakt test (eller laveste positive u-verdi)
        
        # regner ut forkastningssannsynlighet for uvals[uid]
        # og itererer deretter utover i halen til vi finner 
        # approksimert kritiske verdi
        
        tmpSPA = singlesaddle2_BN(uvals[uid],gtilde,muhat,ran)
        while((tmpSPA > alpha) & (uid < length(uvals_true))){
            uid = uid+1
            if(uid == length(uvals_true)){
                # Obs: gjør en korreksjon dersom vi kommer helt ut i halen
                tmpSPA = singlesaddle2_BN(uvals[uid]-0.5,gtilde,muhat,ran)
            }else{
                tmpSPA = singlesaddle2_BN(uvals[uid],gtilde,muhat,ran)
            }
        }
        if(tmpSPA > alpha){
            # Helt ytterst i halen ga ikke forkastning
            signSPA[c] = 0
        }else{
            signSPA[c] = sum(pm[uid:length(pm)])
        }
        critSPA = uid
        
        
        #################
        # SPAcc
        #################
        uid = max(crit-3,which(uvals>0)[1])
        
        tmpSPAcc = singlesaddle2_SCC_BN(uvals[uid],gtilde,muhat,ran)
        while((tmpSPAcc > alpha) & (uid < length(uvals_true))){
            uid = uid+1
            if(uid == length(uvals_true)){
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[uid]-0.5,gtilde,muhat,ran)
            }else{
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[uid],gtilde,muhat,ran)
            }
        }
        
        if(tmpSPAcc > alpha){
            signSPAcc[c] = 0
        }else{
            signSPAcc[c] = sum(pm[uid:length(pm)])
        }
        critSPAcc = uid
        
        
        #################
        # DSPAcc
        #################
        uid = max(crit-3,which(uvals>0)[1])
        
        tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uid],x,muhat,ran)
        while((tmpDSPAcc > alpha) & (uid < length(uvals_true))){
            uid = uid + 1
            if(uid == length(uvals_true)){
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uid]-0.5,x,muhat,ran)
            }else{
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uid],x,muhat,ran)
            }
        }
        
        if(tmpDSPAcc > alpha){
            signDSPAcc[c] = 0
        }else{
            signDSPAcc[c] = sum(pm[uid:length(pm)])
        }
        critDSPAcc = uid
        
        
        #################
        # midPexact
        #################
        uid = max(crit-3,which(uvals>0)[1])
        
        tmpMidP = 0.5*sum(pm[uid:length(pm)]) + 0.5*sum(pm[(uid+1):length(pm)])
        while((tmpMidP > alpha) & (uid < length(uvals))){
            uid = uid + 1
            if(uid == length(uvals)){
                tmpMidP = sum(pm[uid:length(pm)])
            }else{
                tmpMidP = 0.5*sum(pm[uid:length(pm)]) + 0.5*sum(pm[(uid+1):length(pm)])
            }
        }
        
        if(tmpMidP > alpha){
            signMidP[c] = 0
        }else{
            signMidP[c] = sum(pm[uid:length(pm)])
        }
        critMidP = uid
        
        
        
    ############################################################    
    }else if(!is.na(crit)){
        # hvis TRUE:
        # Eksakt test har (i hvert fall) forkasningsregel fra høyre hale
        # men enten to-sidig forkastning eller en-sidig som ligger
        # så "nærme" to-sidig at vi bør sjekke for to-sidige tester
        # med approksimasjonsmetodene
    ############################################################    
        
        #################
        # SPA
        #################
        
        # Først søker vi utover i høyre hale
        # uposid er der vi starter søket
        uposid = max(crit-3,which(uvals>0)[1])
        
        tmpSPA = singlesaddle2_BN(uvals[uposid],gtilde,muhat,ran)
        while(tmpSPA > alpha & uposid < length(uvals)){
            uposid = uposid+1
            if(uposid == length(uvals)){
                tmpSPA = singlesaddle2_BN(uvals[uposid]-0.5,gtilde,muhat,ran)
            }else{
                tmpSPA = singlesaddle2_BN(uvals[uposid],gtilde,muhat,ran)
            }
        }
        
        if(tmpSPA > alpha){
            # Fant ikke noe forkastning i høyre halen
            critSPApos = NA
        }else{
            # Forkastningsregel høyre hale
            critSPApos = uposid
        }
        
        # Deretter søker vi utover i i venstre hale
        # unegid er der vi starter søket
        # vi starter helt ute i halen og ser om det er 
        # forkastning der, eller innover mot null
        unegid = 1
        tmpSPA = singlesaddle2_BN(uvals[unegid]+0.5,gtilde,muhat,ran)
        if(tmpSPA < alpha){ 
        # hvis det går an å forkaste helt i venstre hale, søk innover
            while(tmpSPA < alpha){
                unegid = unegid+1
                tmpSPA = singlesaddle2_BN(uvals[unegid],gtilde,muhat,ran)
            }
            unegid = unegid - 1
            critSPAneg = unegid
            
            if(is.na(critSPApos)){
                # Fant ikke noe forkastning i høyre halen
                # sanns. for type-1 feil bare fra venstre hale
                signSPA[c] = sum(pm[1:critSPAneg])
            }else{
                # sanns. for type-1 feil fra to-sidig forkastningsområde
                signSPA[c] = sum(pm[1:critSPAneg]) + sum(pm[critSPApos:length(pm)])
            }
            
        }else{
        # ikke mulig å forkaste helt i venstre hale, høyre-sidig test
            
            critSPAneg = NA
            if(is.na(critSPApos)){
                # Fant ikke noe forkastning i den høyre halen heller
                signSPA[c] = 0
            }else{
                # sanns. for type-1 feil bare fra høyre hale
                signSPA[c] = sum(pm[critSPApos:length(pm)])
            }
        }
        

        #################
        # SPAcc
        #################
        
        # Fra høyre hale
        uposid = max(crit-3,which(uvals>0)[1])
        tmpSPAcc = singlesaddle2_SCC_BN(uvals[uposid],gtilde,muhat,ran)
        
        while(tmpSPAcc > alpha & uposid < length(uvals)){
            uposid = uposid+1
            if(uposid == length(uvals)){
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[uposid]-0.5,gtilde,muhat,ran)
            }else{
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[uposid],gtilde,muhat,ran)
            }
        }
        
        if(tmpSPAcc > alpha){
            critSPAccpos = NA
        }else{
            critSPAccpos = uposid
        }
        
        # Fra venstre hale
        unegid = 1
        tmpSPAcc = singlesaddle2_SCC_BN(uvals[unegid]+0.5,gtilde,muhat,ran)
        if(tmpSPAcc < alpha){
            while(tmpSPAcc < alpha){
                unegid = unegid+1
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[unegid],gtilde,muhat,ran)
            }
            unegid = unegid-1
            critSPAccneg = unegid
            
            if(is.na(critSPAccpos)){
                signSPAcc[c] = sum(pm[1:critSPAccneg])
            }else{
                signSPAcc[c] = sum(pm[1:critSPAccneg]) + sum(pm[critSPAccpos:length(pm)])
            }
            
        }else{
            critSPAccneg = NA
            if(is.na(critSPAccpos)){
                signSPAcc[c] = 0
            }else{
                signSPAcc[c] = sum(pm[critSPAccpos:length(pm)])
            }
        }
        
        
        #################
        # DSPAcc
        #################
        
        # Fra høyre hale
        uposid = max(crit-3,which(uvals>0)[1])
        tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uposid],x,muhat,ran)
        
        while(tmpDSPAcc > alpha & uposid < length(uvals)){
            uposid = uposid+1
            if(uposid == length(uvals)){
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uposid]-0.5,x,muhat,ran)
            }else{
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uposid],x,muhat,ran)
            }
        }
        if(tmpDSPAcc > alpha){
            critDSPAccpos = NA
        }else{
            critDSPAccpos = uposid
        }
        
        # Fra venstre hale
        unegid = 1
        tmpDSPAcc = doublesaddle2_SCC_BN(uvals[unegid]+0.5,x,muhat,ran)
        if(tmpDSPAcc < alpha){
            while(tmpDSPAcc < alpha){
                unegid = unegid+1
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[unegid],x,muhat,ran)
            }
            unegid = unegid-1
            critDSPAccneg = unegid
            
            if(is.na(critDSPAccpos)){
                signDSPAcc[c] = sum(pm[1:critDSPAccneg])
            }else{
                signDSPAcc[c] = sum(pm[1:critDSPAccneg]) + sum(pm[critDSPAccpos:length(pm)])
            }
            
        }else{
            critDSPAccneg = NA
            if(is.na(critDSPAccpos)){
                signDSPAcc[c] = 0
            }else{
                signDSPAcc[c] = sum(pm[critDSPAccpos:length(pm)])
            }
        }
        
        
        
        #################
        # MidPExact
        #################
        
        # Fra høyre hale
        uposid = max(crit-3,which(uvals>0)[1])
        
        tmpMidP = midp(uposid,uvals,ran,pm)
        
        while(tmpMidP > alpha & uposid < length(uvals)){
            uposid = uposid+1
            tmpMidP = midp(uposid,uvals,ran,pm)
        }
        if(tmpMidP > alpha){
            critMidPpos = NA
        }else{
            critMidPpos = uposid
        }
        
        # Fra venstre hale
        unegid = 1
        tmpMidP = midp(unegid,uvals,ran,pm)
        if(tmpMidP < alpha){
            while(tmpMidP < alpha){
                unegid = unegid+1
                tmpMidP = midp(unegid,uvals,ran,pm)
            }
            unegid = unegid-1
            critMidPneg = unegid
            
            if(is.na(critMidPpos)){
                signMidP[c] = sum(pm[1:critMidPneg])
            }else{
                signMidP[c] = sum(pm[1:critMidPneg]) + sum(pm[critMidPpos:length(pm)])
            }
            
        }else{
            critMidPneg = NA
            if(is.na(critMidPpos)){
                signMidP[c] = 0
            }else{
                signMidP[c] = sum(pm[critMidPpos:length(pm)])
            }
        }
        
        
        
        
    ############################################################ 
    }else{
        # Fordelingen er left-skew og forkastning bare i venstre
        # hale
    ############################################################ 
        crit = t1[3] # eksakt test kritisk verdi i venstre hale
        
        
        #################
        # SPA
        #################
        
        # Starter søket i uid tre hakk inn fra eksakt kritisk verdi
        # eller minste negative u-verdi
        # Søker deretter utover i halen
        uid = min(crit+3,which(uvals<0)[length(which(uvals<0))])
        
        tmpSPA = singlesaddle2_BN(uvals[uid],gtilde,muhat,ran)
        while((tmpSPA > alpha) & (uid > 1)){
            uid = uid-1
            if(uid == 1){
                tmpSPA = singlesaddle2_BN(uvals[uid]+0.5,gtilde,muhat,ran)
            }else{
                tmpSPA = singlesaddle2_BN(uvals[uid],gtilde,muhat,ran)
            }
        }
        
        if (tmpSPA > alpha){
            # fant ikke noen u-verdi som gir forkastning 
            signSPA[c] = 0
        }else{
            signSPA[c] = sum(pm[1:uid])
        }
        critSPA = uid
        
        
        #################
        # SPAcc
        #################
        uid = min(crit+3,which(uvals<0)[length(which(uvals<0))])
        
        tmpSPAcc = singlesaddle2_SCC_BN(uvals[uid],gtilde,muhat,ran)
        while(tmpSPAcc > alpha & uid > 1){
            uid = uid-1
            if(uid == 1){
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[uid]+0.5,gtilde,muhat,ran)
            }else{
                tmpSPAcc = singlesaddle2_SCC_BN(uvals[uid],gtilde,muhat,ran)
            }
        }
        
        if (tmpSPAcc > alpha){
            signSPAcc[c] = 0
        }else{
            signSPAcc[c] = sum(pm[1:uid])
        }
        critSPAcc = uid
        
        
        #################
        # DSPAcc
        #################
        uid = min(crit+3,which(uvals<0)[length(which(uvals<0))])
        
        tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uid],x,muhat,ran)
        while(tmpDSPAcc > alpha & uid > 1){
            uid = uid - 1
            if(uid == 1){
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uid]+0.5,x,muhat,ran)
            }else{
                tmpDSPAcc = doublesaddle2_SCC_BN(uvals[uid],x,muhat,ran)
            }
        }
        
        if (tmpDSPAcc > alpha){
            signSSPAcc[c] = 0
        }else{
            signDSPAcc[c] = sum(pm[1:uid])
        }
        critDSPAcc = uid
        
        
        
        #################
        # MidPExact
        #################
        uid = min(crit+3,which(uvals<0)[length(which(uvals<0))])
        
        tmpMidP = 0.5*sum(pm[1:uid]) + 0.5*sum(pm[1:(uid-1)])
        while(tmpMidP > alpha & uid > 1){
            uid = uid - 1
            if(uid == 1){
                tmpMidP = pm[1]
            }else{
                tmpMidP = 0.5*sum(pm[1:uid]) + 0.5*sum(pm[1:(uid-1)])
            }
        }
        
        if (tmpMidP > alpha){
            signMidP[c] = 0
        }else{
            signMidP[c] = sum(pm[1:uid])
        }
        critMidP = uid
        
    }
    
    
    ############################################################ 
    # normaltilnærming
    ############################################################ 
    sdu = sqrt(mu*(1-mu)*sum(gtilde^2)) # standardavvik
    steps = round(min(uvals[uvals>0]),2) # "desimalene" etter positive heltall
    steps_neg = abs(steps-1) # "desimalene" etter negative heltall
    tmp = qnorm(0.5*alpha,0,sdu,lower.tail = FALSE) # kritisk verdi
    
    # Nå må vi finne tilsvarende u-verdier slik at vi kan regne
    # eksakt sannsynlighet
    
    # høyre hale
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
    
    # venstre hale
    #crit_norm_inv = crit_norm-sign(crit_norm)*ceiling(2*abs(crit_norm))
    #Vi m� benytte oss av forkastningsomr�det til normalfordelingen, se over:
    
    
    if(crit_norm_inv>=ran[1] & crit_norm <= ran[2]){
        indpos = which(round(uvals,2) == round(crit_norm,2))
        indneg = which(round(uvals,2) == round(crit_norm_inv,2))
        signNorm[c] = sum(pm[1:indneg]) + sum(pm[indpos:length(pm)])
    }else if (crit_norm_inv<ran[1] & crit_norm <= ran[2]){
        indpos = which(round(uvals,2) == round(crit_norm,2))
        signNorm[c] = sum(pm[indpos:length(pm)])
    }else if (crit_norm_inv>=ran[1] & crit_norm > ran[2]){
        indneg = which(round(uvals,2) == round(crit_norm_inv,2))
        signNorm[c] = sum(pm[1:indneg])
    }else{
        signNorm[c] = 0
    }
}


signAll = cbind(ncases,
                signExact, 
                signSPA, 
                signSPAcc, 
                signDSPAcc, 
                signNorm,
                signMidP)


nmin = signAll[1,1]
nmax = signAll[dim(signAll)[1],1]

mumin = 0.05
mumax = 0.95
mus = round(seq(mumin,mumax,0.01),2)

# Type-1 feil for hver metode for ulike verdier av mu (under H0)
t1error = matrix(NA,ncol = 6, nrow = length(mus))

# Sannsynlighet for å få et datasett der tester
# ikke holder nivået sitt
probab_invalid = matrix(NA,ncol = 6, nrow = length(mus))

muno = 0

for(mu in mus){
    muno = muno + 1
    yn1 = nmin:nmax # antall kasus vi skal se på
    pn1 = dbinom(yn1,n,mu) # sannsynlighet for antall kasus
    t1error[muno,] = c(sum(signAll[,2]*pn1),
                       sum(signAll[,3]*pn1),
                       sum(signAll[,4]*pn1),
                       sum(signAll[,5]*pn1),
                       sum(signAll[,6]*pn1),
                       sum(signAll[,7]*pn1))
    
    probab_invalid[muno,] = c(
        sum((signAll[,2]>alpha)*pn1),
        sum((signAll[,3]>alpha)*pn1),
        sum((signAll[,4]>alpha)*pn1),
        sum((signAll[,5]>alpha)*pn1),
        sum((signAll[,6]>alpha)*pn1),
        sum((signAll[,7]>alpha)*pn1))
}

minProb = apply(probab_invalid,2,min)
minProbID = apply(probab_invalid,2,which.min)
maxProb = apply(probab_invalid,2,max)
maxProbID = apply(probab_invalid,2,which.max)

prop_invalid = c(sum((signAll[,2]>alpha))/dim(signAll)[1],
                 sum((signAll[,3]>alpha))/dim(signAll)[1], 
                 sum((signAll[,4]>alpha))/dim(signAll)[1],
                 sum((signAll[,5]>alpha))/dim(signAll)[1],
                 sum((signAll[,6]>alpha))/dim(signAll)[1],
                 sum((signAll[,7]>alpha))/dim(signAll)[1])


which(signAll[,4]>alpha)
signAll[c(401,589),1]



colnames(t1error) = c("Exact","SPA","SPACC","DSPACC","Norm","MidP")
colnames(probab_invalid) = c("Exact","SPA","SPACC","DSPACC","Norm","MidP")

t1error = as.data.frame(t1error)
probab_invalid = as.data.frame(probab_invalid)


#df = data.frame(t1error)
#df = melt(df)
df= data.frame(mu = rep(mus,5))
df$t1errors = c(t1error$Exact,t1error$SPA,t1error$SPACC,t1error$DSPACC,t1error$Norm)
df$Method = c(rep("Exact",nrow(t1error)),rep("ESPA",nrow(t1error)),rep("ESPA-CC",nrow(t1error)),rep("DSPA-CC",nrow(t1error)),rep("Norm",nrow(t1error)))


df2= data.frame(mu = rep(mus,5))
df2$probab_invalids = c(probab_invalid$Exact,probab_invalid$SPA,probab_invalid$SPACC,probab_invalid$DSPACC,probab_invalid$Norm)
df2$Method = c(rep("Exact",nrow(t1error)),rep("ESPA",nrow(t1error)),rep("ESPA-CC",nrow(t1error)),rep("DSPA-CC",nrow(t1error)),rep("Norm",nrow(t1error)))

if(r == 1){
p1 = ggplot(df, aes(x = mu, y = t1errors, color = Method)) +
    theme_classic()+
    geom_line() +
    geom_abline(intercept = alpha,slope = 0) +
    scale_y_continuous("Overall type I error probabilities") +
    coord_cartesian(ylim=c(0, 0.063)) +
    scale_color_viridis(option="turbo",discrete = TRUE) +
    ggtitle(paste0("Signif. level ",alpha)) +
    labs(x=expression(mu))+
    scale_x_continuous(breaks = seq(0,1,0.1))

p1

#df3 = data.frame(variable = unique(df2$variable),
#                 prop_invalid = prop_invalid)

p2 = ggplot(df2, aes(x = mu, y = probab_invalids,
                     color = Method)) +
    theme_classic()+
    geom_line() +
    coord_cartesian(ylim= c(0,1)) +
    scale_y_continuous("Probability of conditionally invalid test") +
    scale_color_viridis(option="turbo",discrete = TRUE) +
    ggtitle(paste0("Signif. level ",alpha)) +
    labs(x=expression(mu))+
    scale_x_continuous(breaks = seq(0.1,0.9,0.1)) #+ 
    #geom_hline(data = df3, 
    #           aes(yintercept = prop_invalid, color = variable),
    #           lty =2)



p2
}else{
    p3 = ggplot(df, aes(x = mu, y = t1errors, color = Method)) +
        theme_classic()+
        geom_line() +
        geom_abline(intercept = alpha,slope = 0) +
        scale_y_continuous("Overall type I error probabilities") +
        coord_cartesian(ylim=c(0, 1e-4)) +
        scale_color_viridis(option="turbo",discrete = TRUE) +
        ggtitle(paste0("Signif. level ",alpha)) +
        labs(x=expression(mu))+
        scale_x_continuous(breaks = seq(0,1,0.1))
    
    p3
    
    #df3 = data.frame(variable = unique(df2$variable),
    #                 prop_invalid = prop_invalid)
    
    p4 = ggplot(df2, aes(x = mu, y = probab_invalids,
                         color = Method)) +
        theme_classic()+
        geom_line() +
        coord_cartesian(ylim= c(0,1)) +
        scale_y_continuous("Probability of conditionally invalid test") +
        scale_color_viridis(option="turbo",discrete = TRUE) +
        ggtitle(paste0("Signif. level ",alpha)) +
        labs(x=expression(mu))+
        scale_x_continuous(breaks = seq(0.1,0.9,0.1)) #+ 
    #geom_hline(data = df3, 
    #           aes(yintercept = prop_invalid, color = variable),
    #           lty =2)
    
    
    
    p4    
    
}

}

ggarrange(p1,p3,p2,p4,byrow=FALSE,nrow = 2,ncol = 2,common.legend = TRUE)



# plot(x = signAll[,1], y = signAll[,3], type = "l")
# lines(x = signAll[,1], y = signAll[,6], type = "l", col = "red")
# abline(h = alpha)
# 
# plot(x = signAll[,1], y = signAll[,6], type = "l")
# lines(x = signAll[,1], y = signAll[,3], type = "l", col = "red")
# abline(h = alpha)
# 
# 
# 
# minProb = c()
# maxProb = c()





