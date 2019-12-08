#input parameters
# 3.1 a Transition probabilities for inflectra
p_ACR20        <- (72.6-39.5)/100      # probability of ACR50>x>ACR20 at the beginning of a cycle for each 
p_ACR50        <- (39.5-16.5)/100      # probability of ACR70>x>ACR50 at the beginning of a cycle for each 
p_ACR70        <- 16.5/100     # probability of x>ACR70 at the beginning of a cycle for each 
p_nonresp <- (100-72.6)/100  # probability if x<ACR20
p_discontinue <- (8/302)*(7.5/40) #probability of discontinuation due to an adverse event - model cost accordingly at wk 30/ 7 months 28 out of 302 discontinued due to ae

# 3.1 a Transition probabilities for remicade
#specify difference
ACR20diffL<--.01
ACR20diffH<-.15
ACR50diffL<--.03
ACR50diffH<-.14
ACR70diffL<--.03
ACR70diffH<-.09

ACR20diff<-.07
ACR50diff<-.05
ACR70diff<-.03

p_ACR20r        <- (65.3-33.9)/100      # probability of ACR50>x>ACR20 at the beginning of a cycle for each 
p_ACR50r        <- (33.9-13.5)/100      # probability of ACR70>x>ACR50 at the beginning of a cycle for each 
p_ACR70r        <- 13.5/100     # probability of x>ACR70 at the beginning of a cycle for each 
p_nonrespr <- (100-65.3)/100  # probability if x<ACR20
p_discontinuer <- (22/304)*(7.5/40) #probability of discontinuation due to an adverse event - model cost accordingly at wk 30/ 7 months 28 out of 302 discontinued due to ae


#Define functions for sensitivity analysis
Probs3 <- function(M_it, v_haq, v_age_at,t,ACR20diff_at=ACR20diff,ACR50diff_at=ACR50diff,ACR70diff_at=ACR70diff) { 
  # Arguments:
  # M_it: health state occupied by individual i at cycle t (character variable)
  # v_haq: vector with haq score at cycle t
  # v_age: vector with age at cycle t
  # t:     current cycle 
  # Returns: 
  #   transition probabilities for that cycle
  v_age <- data.frame(agecol=v_age_at)
  
  m_p_it           <- matrix(0, nrow = n_s+1, ncol = n_i)  # create matrix of state transition probabilities
  rownames(m_p_it) <- c("baseline",v_n)# give the state names to the rows
  
  # lookup baseline probability and rate of dying based on age and haq
  p_death_age <- inner_join(v_age, m_p_HD, by=c("agecol") ) #lookup prob death based on age 
  p_death_age_haq <- p_death_age
  p_death_age_haq$mortality<-p_death_age_haq$mortality*1.33^(v_haq)
  p_mort<-p_death_age_haq$mortality
  
  p_ACR20new<-p_ACR20-ACR20diff_at
  p_ACR50new<-p_ACR50-ACR50diff_at
  p_ACR70new<-p_ACR70-ACR70diff_at
  p_nonrespNew<-1-(p_ACR20new+p_ACR50new+p_ACR70new)
  

  
  for (i in 1:n_i){
    m_p_it[1,i]<-0 #no one moves to baseline
    m_p_it[2,i]<-ifelse(M_it[i]=="baseline",(1-(1-p_mort[i])*p_discontinue)*p_ACR20new,m_p_it[2,i])#moving to ACR20 from baseline
    m_p_it[3,i]<-ifelse(M_it[i]=="baseline",(1-(1-p_mort[i])*p_discontinue)*p_ACR50new,m_p_it[3,i])#moving to ACR50 from baseline
    m_p_it[4,i]<-ifelse(M_it[i]=="baseline",(1-(1-p_mort[i])*p_discontinue)*p_ACR70new,m_p_it[4,i])#moving to ACR50 from baseline
    m_p_it[5,i]<-ifelse(M_it[i]=="baseline",(1-(1-p_mort[i])*p_discontinue)*p_nonrespNew,m_p_it[5,i])#moving to nonresp from baseline
    m_p_it[6,i]<-ifelse(M_it[i]=="baseline",(1-p_mort[i])*p_discontinue,m_p_it[6,i])#moving to discontinue from baseline
    m_p_it[7,i]<-ifelse(M_it[i]=="baseline",0,m_p_it[7,i])#moving to lastline from baseline
    m_p_it[8,i]<-ifelse(M_it[i]=="baseline",p_mort[i],m_p_it[8,i])#moving to death from baseline
    
    m_p_it[2,i]<-ifelse(M_it[i]=="ACR20",(1-(1-p_mort[i])*p_discontinue)*p_ACR20new,m_p_it[2,i])
    m_p_it[3,i]<-ifelse(M_it[i]=="ACR20",(1-(1-p_mort[i])*p_discontinue)*p_ACR50new,m_p_it[3,i])
    m_p_it[4,i]<-ifelse(M_it[i]=="ACR20",(1-(1-p_mort[i])*p_discontinue)*p_ACR70new,m_p_it[4,i])
    m_p_it[5,i]<-ifelse(M_it[i]=="ACR20",(1-(1-p_mort[i])*p_discontinue)*p_nonrespNew,m_p_it[5,i])#moving to nonresp from baseline
    m_p_it[6,i]<-ifelse(M_it[i]=="ACR20",(1-p_mort[i])*p_discontinue,m_p_it[6,i])#moving to discontinue from baseline
    m_p_it[7,i]<-ifelse(M_it[i]=="ACR20",0,m_p_it[7,i])#moving to lastline from baseline
    m_p_it[8,i]<-ifelse(M_it[i]=="ACR20",p_mort[i],m_p_it[8,i])#moving to death from baseline
    
    m_p_it[2,i]<-ifelse(M_it[i]=="ACR50",(1-(1-p_mort[i])*p_discontinue)*p_ACR20new,m_p_it[2,i])
    m_p_it[3,i]<-ifelse(M_it[i]=="ACR50",(1-(1-p_mort[i])*p_discontinue)*p_ACR50new,m_p_it[3,i])
    m_p_it[4,i]<-ifelse(M_it[i]=="ACR50",(1-(1-p_mort[i])*p_discontinue)*p_ACR70new,m_p_it[4,i])
    m_p_it[5,i]<-ifelse(M_it[i]=="ACR50",(1-(1-p_mort[i])*p_discontinue)*p_nonrespNew,m_p_it[5,i])#moving to nonresp from baseline
    m_p_it[6,i]<-ifelse(M_it[i]=="ACR50",(1-p_mort[i])*p_discontinue,m_p_it[6,i])#moving to discontinue from baseline
    m_p_it[7,i]<-ifelse(M_it[i]=="ACR50",0,m_p_it[7,i])#moving to lastline from baseline
    m_p_it[8,i]<-ifelse(M_it[i]=="ACR50",p_mort[i],m_p_it[8,i])#moving to death from baseline
    
    m_p_it[2,i]<-ifelse(M_it[i]=="ACR70",(1-(1-p_mort[i])*p_discontinue)*p_ACR20new,m_p_it[2,i])
    m_p_it[3,i]<-ifelse(M_it[i]=="ACR70",(1-(1-p_mort[i])*p_discontinue)*p_ACR50new,m_p_it[3,i])
    m_p_it[4,i]<-ifelse(M_it[i]=="ACR70",(1-(1-p_mort[i])*p_discontinue)*p_ACR70new,m_p_it[4,i])
    m_p_it[5,i]<-ifelse(M_it[i]=="ACR70",(1-(1-p_mort[i])*p_discontinue)*p_nonrespNew,m_p_it[5,i])#moving to nonresp from baseline
    m_p_it[6,i]<-ifelse(M_it[i]=="ACR70",(1-p_mort[i])*p_discontinue,m_p_it[6,i])#moving to discontinue from baseline
    m_p_it[7,i]<-ifelse(M_it[i]=="ACR70",0,m_p_it[7,i])#moving to lastline from baseline
    m_p_it[8,i]<-ifelse(M_it[i]=="ACR70",p_mort[i],m_p_it[8,i])#moving to death from baseline
    
    m_p_it[2,i]<-ifelse(M_it[i]=="nonresp",0 ,m_p_it[2,i])#moving to all other states from nonresp
    m_p_it[3,i]<-ifelse(M_it[i]=="nonresp",0 ,m_p_it[3,i])#moving to all other states from nonresp
    m_p_it[4,i]<-ifelse(M_it[i]=="nonresp",0 ,m_p_it[4,i])#moving to all other states from nonresp
    m_p_it[5,i]<-ifelse(M_it[i]=="nonresp",0 ,m_p_it[5,i])#moving to all other states from nonresp
    m_p_it[6,i]<-ifelse(M_it[i]=="nonresp",0 ,m_p_it[6,i])#moving to all other states from nonresp
    m_p_it[7,i]<-ifelse(M_it[i]=="nonresp",1-p_mort[i],m_p_it[7,i])#moving to lastline from nonresp
    m_p_it[8,i]<-ifelse(M_it[i]=="nonresp",p_mort[i],m_p_it[8,i])#moving to death from nonresp
    
    m_p_it[2,i]<-ifelse(M_it[i]=="discontinue",0 ,m_p_it[2,i])#moving to all other states from nonresp
    m_p_it[3,i]<-ifelse(M_it[i]=="discontinue",0 ,m_p_it[3,i])#moving to all other states from nonresp
    m_p_it[4,i]<-ifelse(M_it[i]=="discontinue",0 ,m_p_it[4,i])#moving to all other states from nonresp
    m_p_it[5,i]<-ifelse(M_it[i]=="discontinue",0 ,m_p_it[5,i])#moving to all other states from nonresp
    m_p_it[6,i]<-ifelse(M_it[i]=="discontinue",0 ,m_p_it[6,i])#moving to all other states from nonresp
    m_p_it[7,i]<-ifelse(M_it[i]=="discontinue", 1-p_mort[i],m_p_it[7,i])#moving to lastline from nonresp
    m_p_it[8,i]<-ifelse(M_it[i]=="discontinue", p_mort[i],m_p_it[8,i])#moving to death from nonresp
    
    m_p_it[2,i]<-ifelse(M_it[i]=="last_line",0 ,m_p_it[2,i])#moving to all other states from lastnine
    m_p_it[3,i]<-ifelse(M_it[i]=="last_line",0 ,m_p_it[3,i])#moving to all other states from lastnine
    m_p_it[4,i]<-ifelse(M_it[i]=="last_line",0 ,m_p_it[4,i])#moving to all other states from lastnine
    m_p_it[5,i]<-ifelse(M_it[i]=="last_line",0 ,m_p_it[5,i])#moving to all other states from lastnine
    m_p_it[6,i]<-ifelse(M_it[i]=="last_line",0 ,m_p_it[6,i])#moving to all other states from lastnine
    m_p_it[7,i]<-ifelse(M_it[i]=="last_line",1-p_mort[i] ,m_p_it[7,i])#moving to all other states from lastnine
    m_p_it[8,i]<-ifelse(M_it[i]=="last_line",p_mort[i],m_p_it[8,i])#moving to death from nonresp
    
    m_p_it[2,i]<-ifelse(M_it[i]=="Dead",0 ,m_p_it[2,i])#moving to all other states from lastnine
    m_p_it[3,i]<-ifelse(M_it[i]=="Dead",0 ,m_p_it[3,i])#moving to all other states from lastnine
    m_p_it[4,i]<-ifelse(M_it[i]=="Dead",0 ,m_p_it[4,i])#moving to all other states from lastnine
    m_p_it[5,i]<-ifelse(M_it[i]=="Dead",0 ,m_p_it[5,i])#moving to all other states from lastnine
    m_p_it[6,i]<-ifelse(M_it[i]=="Dead",0 ,m_p_it[6,i])#moving to all other states from lastnine
    m_p_it[7,i]<-ifelse(M_it[i]=="Dead",0 ,m_p_it[7,i])#moving to all other states from lastnine
    m_p_it[8,i]<-ifelse(M_it[i]=="Dead",1,m_p_it[8,i])#moving to death from nonresp
    m_p_it= apply(m_p_it,2,function(x){x/sum(x)})
    
  }
  return(t(m_p_it))
}
MicroSim3 <- function(n_i, seed = 1,ACR20diffSim=ACR20diff,ACR50diffSim=ACR50diff,ACR70diffSim=ACR70diff) {
  # Arguments:  
  # n_i:     number of individuals
  # df_X     data frame with individual data 
  # seed:    default is 1
  
  set.seed(seed) # set the seed
  
  n_s <- length(v_n) # the number of health states
  
  # create three matrices called m_M, m_C and m_E
  # number of rows is equal to the n_i, the number of columns is equal to n_t (the initial state and all the n_t cycles)
  # m_M is used to store the health state information over time for every individual
  # m_C is used to store the costs information over time for every individual
  # m_E is used to store the effects information over time for every individual
  # m_age is used to store the age information for each individual
  m_M <- m_C <- m_E <-  matrix(nrow = n_i, ncol = n_t + 1, 
                               dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                               paste("cycle", 0:n_t, sep = " ")))  
  
  m_M[, 1] <- v_M_init         # initial health state for individual i
  v_Ts     <- v_Ts_init        # initialize time since illness onset for individual i
  v_Ts_cont     <- v_Ts_init
  
  v_age_sim<-v_age_init #initialize age
  v_age_cont<-v_age_init#initialize cont_age
  
  v_last_dur<-v_last_dur_init
  v_last_dur_cont<-v_last_dur_init

  
  v_HAQ_at<-v_HAQ_init
  m_C[, 1] <- Costs2(M_it=m_M[, 1])   # costs accrued individual i during cycle 0
  m_E[, 1] <- matrix(Effs2(M_it=m_M[, 1],HAQ_at=v_HAQ_init,v_TS_at=v_Ts_init,v_last_dur_at=v_last_dur_init,age_at=v_age_sim))    # QALYs accrued individual i during cycle 0
  # open a loop for time running cycles 1 to n_t 
  for (t in 1:n_t) {
    v_p <- Probs3(M_it=m_M[, t], v_haq=v_HAQ_at, v_age=v_age_sim, t,ACR20diff_at=ACR20diffSim,ACR50diff_at=ACR50diffSim,ACR70diff_at=ACR70diffSim)  # calculate the transition probabilities for the cycle based on health state t
    m_M[, t + 1]  <- samplev(v_p, 1)       # sample the current health state and store that state in matrix m_M 
    m_C[, t + 1]  <- Costs2(M_it=m_M[, t + 1])   # calculate costs per individual during cycle t + 1
    m_E[, t + 1]  <- matrix(Effs2 (M_it=m_M[, t + 1],v_TS_at = v_Ts,HAQ_at =v_HAQ_at,v_last_dur_at = v_last_dur,age_at=v_age_sim))   # calculate QALYs per individual during cycle t + 1
    
    ifelse(v_age_cont<110,v_age_cont+.144,v_age_cont)
    v_age_Sim<-round(v_age_cont)
    #update HAQ
    v_HAQ_at<-if_else(m_M[, t + 1] == "ACR70" & memory70!=1,v_HAQ_at-1.07,v_HAQ_at)
    v_HAQ_at<-if_else(m_M[, t + 1] == "ACR50"& memory50!=1,v_HAQ_at-.76,v_HAQ_at)
    v_HAQ_at<-if_else(m_M[, t + 1] == "ACR20"& memory20!=1,v_HAQ_at-.44,v_HAQ_at)
    v_HAQ_at<-if_else(m_M[, t + 1] == "nonresp"& memorynonresp!=1,v_HAQ_at-.11,v_HAQ_at)
    v_HAQ_at<-if_else(m_M[, t + 1] == "last_line"&v_last_dur==0,v_HAQ_at+.11,v_HAQ_at)
    
    memory70<-if_else(m_M[, t + 1] == "ACR70",1,0)
    memory50<-if_else(m_M[, t + 1] == "ACR50",1,0)
    memory20<-if_else(m_M[, t + 1] == "ACR20",1,0)
    memorynonresp<-if_else(m_M[, t + 1] == "nonresp",1,0)
    
    v_Ts_cont<-if_else(m_M[, t + 1] != "Dead", v_Ts_cont + .144, v_Ts_cont) # update time since illness onset for t + 1 
    v_last_dur_cont<-if_else(m_M[, t + 1] == "last_line", v_last_dur_cont+.144, v_last_dur_cont)
    v_Ts<-round(v_Ts_cont)
    v_last_dur<-round(v_last_dur_cont)
    
    # Display simulation progress
    if(t/(n_t/10) == round(t/(n_t/10), 0)) { # display progress every 10%
      cat('\r', paste(t/n_t * 100, "% done", sep = " "))
    }
    
  } # close the loop for the time points 
  
  # calculate  
  tc <- m_C%*%(v_dwc[1:61])     #total (discounted) cost per individual
  te <- m_E%*%(v_dwe[1:61])      #total (discounted) QALYs per individual 
  tc_hat <- mean(tc)     # average (discounted) cost 
  te_hat <- mean(te)     # average (discounted) QALYs
  
  # store the results from the simulation in a list
  results <- list(m_M = m_M, m_C = m_C, m_E = m_E, tc = tc , te = te, tc_hat = tc_hat, te_hat = te_hat)   
  return(results)  # return the results
} # end of the MicroSim function  

#run microsimulations
outcomes_Base<-MicroSim3(n_i)
outcomes_sens20L<-MicroSim3(n_i,ACR20diffSim=ACR20diffL)
outcomes_sens20H<-MicroSim3(n_i,ACR20diffSim=ACR20diffH)
outcomes_sens50L<-MicroSim3(n_i,ACR50diffSim=ACR50diffL)
outcomes_sens50H<-MicroSim3(n_i,ACR50diffSim=ACR50diffH)
outcomes_sens70L<-MicroSim3(n_i,ACR70diffSim=ACR70diffL)
outcomes_sens70H<-MicroSim3(n_i,ACR70diffSim=ACR70diffH)
outcomes_inf<-MicroSim(n_i)

#### 08 Cost Effectiveness Analysis ####
# store the mean costs of each strategy in a new variable C (vector of costs)
v_CBase <- c(outcomes_inf$tc_hat, outcomes_Base$tc_hat)
v_C20L <- c(outcomes_inf$tc_hat, outcomes_sens20L$tc_hat)
v_C20H <- c(outcomes_inf$tc_hat, outcomes_sens20H$tc_hat)
v_C50L <- c(outcomes_inf$tc_hat, outcomes_sens50L$tc_hat)
v_C50H <- c(outcomes_inf$tc_hat, outcomes_sens50H$tc_hat)
v_C70L <- c(outcomes_inf$tc_hat, outcomes_sens70L$tc_hat)
v_C70H <- c(outcomes_inf$tc_hat, outcomes_sens70H$tc_hat)

# store the mean QALYs of each strategy in a new variable E (vector of effects)
v_EBase <- c(outcomes_inf$te_hat, outcomes_Base$te_hat)
v_E20L <- c(outcomes_inf$te_hat, outcomes_sens20L$te_hat)
v_E20H <- c(outcomes_inf$te_hat, outcomes_sens20H$te_hat)
v_E50L <- c(outcomes_inf$te_hat, outcomes_sens50L$te_hat)
v_E50H <- c(outcomes_inf$te_hat, outcomes_sens50H$te_hat)
v_E70L <- c(outcomes_inf$te_hat, outcomes_sens70L$te_hat)
v_E70H <- c(outcomes_inf$te_hat, outcomes_sens70H$te_hat)

# strategy names
v_names_str <- c("INF", "REM") 

# use dampack to calculate the ICER
ICERbase<-calculate_icers(cost= v_CBase,effect= v_EBase,strategies = v_names_str)
ICERbase<-(v_CBase[1]-v_CBase[2])/(v_EBase[1]/v_EBase[2])
ICER20L<-calculate_icers(cost= v_C20L,effect= v_E20L,strategies = v_names_str)
ICER20L<-(v_C20L[1]-v_C20L[2])/(v_E20L[1]/v_E20L[2])

ICER20H<-calculate_icers(cost= v_C20H,effect= v_E20H,strategies = v_names_str)
ICER20H<-(v_C20H[1]-v_C20H[2])/(v_E20H[1]/v_E20H[2])

ICER50L<-calculate_icers(cost= v_C50L,effect= v_E50L,strategies = v_names_str)
ICER50L<-(v_C50L[1]-v_C50L[2])/(v_E50L[1]/v_E50L[2])

ICER50H<-calculate_icers(cost= v_C50H,effect= v_E50H,strategies = v_names_str)
ICER50H<-(v_C50H[1]-v_C50H[2])/(v_E50H[1]/v_E50H[2])

ICER70L<-calculate_icers(cost= v_C70L,effect= v_E70L,strategies = v_names_str)
ICER70L<-(v_C70L[1]-v_C70L[2])/(v_E70L[1]/v_E70L[2])

ICER70H<-calculate_icers(cost= v_C70H,effect= v_E70H,strategies = v_names_str)
ICER70H<-(v_C70H[1]-v_C70H[2])/(v_E70H[1]/v_E70H[2])

#plot

plot(density(outcomes_Base$tc), main = paste("REM Total cost per person"), xlab = "Cost ($)")
plot(density(outcomes_Base$te), main = paste("REM Total QALYs per person"), xlab = "QALYs")

