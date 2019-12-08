
###############  RA biologic Microsimulation  #############
# Includes: 
# individual characteristics: previous therapy
#next therpy depends on acr and previous therapy
################################################################################

################################################################################
# Please cite our publications when using this code
# darthworkgroup.com 
## Jalal H, et al. An Overview of R in Health Decision Sciences. 
# Med. Decis. Making. 2017; 37(3): 735-746. 
## Krijkamp EM, et al. Microsimulation modeling for health decision sciences 
# using R: a tutorial. Med. Decis. Making. 2018; 38(3): 400-422.

################################################################################
rm(list = ls())  # Delete everything that is in R's memory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory as the folder where the course material is stored
#### 01 Load packages ####n
library(dplyr)     # load plyr including the useful join function
library(dampack)  # for CEA, PSA, and visualization
library(reshape2) # for data manipulation
library(ggplot2)
#### 02 Load Functions ####
source("Functions.R")

#### 03 Input Model Parameters ####
set.seed(1)  # set the seed  

# 3.1 Model structure
v_n   <- c("ACR20","ACR50","ACR70", "nonresp", "discontinue","last_line","Dead")          # vector with state names
n_s   <- length(v_n)                           # number of states
n_t   <- 60                                    # number of cycles
n_i   <- 1000                                 # number of individuals
d_r   <- 0.03                                  # discount rate of 3% per cycle
v_dwe <- v_dwc <- 1 / ((1 + d_r) ^ (0:n_t))    # discount weight 

memory70<-2
memory50<-2
memory20<-2
memorynonresp<-2
# 3.1 a Transition probabilities for inflectra
p_ACR20        <- (72.6-39.5)/100      # probability of ACR50>x>ACR20 at the beginning of a cycle for each 
p_ACR50        <- (39.5-16.5)/100      # probability of ACR70>x>ACR50 at the beginning of a cycle for each 
p_ACR70        <- 16.5/100     # probability of x>ACR70 at the beginning of a cycle for each 
p_nonresp <- (100-72.6)/100  # probability if x<ACR20
p_discontinue <- (8/302)*(7.5/40) #probability of discontinuation due to an adverse event - model cost accordingly at wk 30/ 7 months 28 out of 302 discontinued due to ae

# 3.1 a Transition probabilities for remicade
p_ACR20r        <- (65.3-33.9)/100      # probability of ACR50>x>ACR20 at the beginning of a cycle for each 
p_ACR50r        <- (33.9-13.5)/100      # probability of ACR70>x>ACR50 at the beginning of a cycle for each 
p_ACR70r        <- 13.5/100     # probability of x>ACR70 at the beginning of a cycle for each 
p_nonrespr <- (100-65.3)/100  # probability if x<ACR20
p_discontinuer <- (22/304)*(7.5/40) #probability of discontinuation due to an adverse event - model cost accordingly at wk 30/ 7 months 28 out of 302 discontinued due to ae

#3.1 b take mortality from life table
lowest_age     <- 40                       # age at baseline
max_age <- 111                       # maximum age of follow up
lt_usa_2005 <- read.csv("HMD_USA_Mx_2015.csv")
v_r_HD <- lt_usa_2005 %>% 
  filter(Age >= lowest_age & Age <= (max_age-1)) %>%
  select(Total) %>%
  as.matrix()
p_HD   <- 1 - exp(- v_r_HD)  # probability of dying from ages 50 to 110 when healthy
m_p_HD<-data.frame(agecol=as.integer(40:110), mortality = p_HD)
colnames(m_p_HD)<-c("agecol","mortality")

#3.2 cost inputs
#3.2 a drug costs
c_inf<-946.28*3/100 #inflectra unit: unit wac 100mg / 3mg/kg weeks 0,2, 6 and every 8 weeks
c_rem<-1167.82*3/100 #remicade unit: unit wac 100mg / weeks 3mg/kg weeks 0,2, 6 and every 8 weeks
c_mtx<-0.02*325.13*7.5 #methotrexate as last line therapy 2.5mg / alone 20mg/wk, coadministration mtx 15.6 (3.1) for inf and 15.6 mg/wk 
c_mtxCo<-0.0156*325.13*7.5 #coadministration mtx 15.6 mg/wk 

#3.2 b admin and monitor cost inputs
c_admin<-136.41+28.64+73.4 #1st and second hour iv and cost per office visit
c_monitor<-(73.4+84.95+4*7.63+4*10.67)*(7.5/52.1429) #1 office 1TB, 4Liver, 4blood count, adjusted for 7.5 wks cycle
#3.3 adverse event inputs
c_infect<-13747
c_tuberculosis<-12220
c_ae_inf<-(8)*c_infect/(8+16)+(16)*c_tuberculosis/(8+16)
c_ae_rem<-(11)*c_infect/(11+14)+(14)*c_tuberculosis/(11+14)
u_tuberculosis<-.156*2/7.5
u_infect<-.156/7.5
u_ae_inf<-(8)*u_infect/(8+16)+(16)*u_tuberculosis/(8+16)
u_ae_rem<-(11)*u_infect/(11+14)+(14)*u_tuberculosis/(11+14)


#### 04 Sample individual level characteristics ####
#### 04.1 Static characteristics ####
v_sex    <- sample(x = c("Female", "Male"), prob = c(0.79, .21), size = n_i, replace = TRUE) # randomly sample the sex of an individual (50% female)
df_x<- data.frame(ID = 1:n_i, Sex = v_sex,weight=0)  # create matrix of state transition probabilities
# lookup sex and update weight taken from https://www.cdc.gov/nchs/data/nhsr/nhsr122-508.pdf
df_x2<-df_x %>%
    rowwise() %>%
    mutate(weight =  round(rnorm(1, mean = ifelse(Sex == "Female", 170.5/2.205 , 197.8/2.205), ifelse(Sex == "Female", 1.7/2.205 , 1.9/2.205)), 1))
df_x3<-df_x2
df_x3$sex_binary[df_x3[,2]=="Male"]<-1
df_x3$sex_binary[df_x3[,2]=="Female"]<-0


#### 04.2 Dynamic characteristics 
# Specify the initial health baseline state of the individuals 
v_M_init  <- rep("baseline", times = n_i)   
v_HAQ_init <- rnorm(n_i,mean=1.6,sd=.6)   # from plantera a vector with the number of previous tmt at the start of the model 
v_last_dur_init<-rep(0,n_i) #duration of last line therapy
v_age_init <- sample(40:60, n_i, T)
v_Ts_init <- rep(0, n_i)  # a vector with the time of being sick at the start of the model 
df_age_init <- data.frame(agecol=v_age_init)


#### 05.1 Probability function ####
# The Probs function that updates the transition probabilities of every cycle is shown below.For tmt or no tmt

Probs <- function(M_it, v_haq, v_age_at,t) { 
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
  
  for (i in 1:n_i){
    m_p_it[1,i]<-0 #no one moves to baseline
    m_p_it[2,i]<-ifelse(M_it[i]=="baseline",(1-(1-p_mort[i])*p_discontinue)*p_ACR20,m_p_it[2,i])#moving to ACR20 from baseline
    m_p_it[3,i]<-ifelse(M_it[i]=="baseline",(1-(1-p_mort[i])*p_discontinue)*p_ACR50,m_p_it[3,i])#moving to ACR50 from baseline
    m_p_it[4,i]<-ifelse(M_it[i]=="baseline",(1-(1-p_mort[i])*p_discontinue)*p_ACR70,m_p_it[4,i])#moving to ACR50 from baseline
    m_p_it[5,i]<-ifelse(M_it[i]=="baseline",(1-(1-p_mort[i])*p_discontinue)*p_nonresp,m_p_it[5,i])#moving to nonresp from baseline
    m_p_it[6,i]<-ifelse(M_it[i]=="baseline",(1-p_mort[i])*p_discontinue,m_p_it[6,i])#moving to discontinue from baseline
    m_p_it[7,i]<-ifelse(M_it[i]=="baseline",0,m_p_it[7,i])#moving to lastline from baseline
    m_p_it[8,i]<-ifelse(M_it[i]=="baseline",p_mort[i],m_p_it[8,i])#moving to death from baseline
    
    m_p_it[2,i]<-ifelse(M_it[i]=="ACR20",(1-(1-p_mort[i])*p_discontinue)*p_ACR20,m_p_it[2,i])
    m_p_it[3,i]<-ifelse(M_it[i]=="ACR20",(1-(1-p_mort[i])*p_discontinue)*p_ACR50,m_p_it[3,i])
    m_p_it[4,i]<-ifelse(M_it[i]=="ACR20",(1-(1-p_mort[i])*p_discontinue)*p_ACR70,m_p_it[4,i])
    m_p_it[5,i]<-ifelse(M_it[i]=="ACR20",(1-(1-p_mort[i])*p_discontinue)*p_nonresp,m_p_it[5,i])#moving to nonresp from baseline
    m_p_it[6,i]<-ifelse(M_it[i]=="ACR20",(1-p_mort[i])*p_discontinue,m_p_it[6,i])#moving to discontinue from baseline
    m_p_it[7,i]<-ifelse(M_it[i]=="ACR20",0,m_p_it[7,i])#moving to lastline from baseline
    m_p_it[8,i]<-ifelse(M_it[i]=="ACR20",p_mort[i],m_p_it[8,i])#moving to death from baseline
    
    m_p_it[2,i]<-ifelse(M_it[i]=="ACR50",(1-(1-p_mort[i])*p_discontinue)*p_ACR20,m_p_it[2,i])
    m_p_it[3,i]<-ifelse(M_it[i]=="ACR50",(1-(1-p_mort[i])*p_discontinue)*p_ACR50,m_p_it[3,i])
    m_p_it[4,i]<-ifelse(M_it[i]=="ACR50",(1-(1-p_mort[i])*p_discontinue)*p_ACR70,m_p_it[4,i])
    m_p_it[5,i]<-ifelse(M_it[i]=="ACR50",(1-(1-p_mort[i])*p_discontinue)*p_nonresp,m_p_it[5,i])#moving to nonresp from baseline
    m_p_it[6,i]<-ifelse(M_it[i]=="ACR50",(1-p_mort[i])*p_discontinue,m_p_it[6,i])#moving to discontinue from baseline
    m_p_it[7,i]<-ifelse(M_it[i]=="ACR50",0,m_p_it[7,i])#moving to lastline from baseline
    m_p_it[8,i]<-ifelse(M_it[i]=="ACR50",p_mort[i],m_p_it[8,i])#moving to death from baseline
    
    m_p_it[2,i]<-ifelse(M_it[i]=="ACR70",(1-(1-p_mort[i])*p_discontinue)*p_ACR20,m_p_it[2,i])
    m_p_it[3,i]<-ifelse(M_it[i]=="ACR70",(1-(1-p_mort[i])*p_discontinue)*p_ACR50,m_p_it[3,i])
    m_p_it[4,i]<-ifelse(M_it[i]=="ACR70",(1-(1-p_mort[i])*p_discontinue)*p_ACR70,m_p_it[4,i])
    m_p_it[5,i]<-ifelse(M_it[i]=="ACR70",(1-(1-p_mort[i])*p_discontinue)*p_nonresp,m_p_it[5,i])#moving to nonresp from baseline
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

testprobs<-Probs(M_it=v_M_test_1,v_haq=v_HAQ_init,v_age=v_age_init,3)
rowSums(testprobs)

#### 05.2 Cost function ####
# The Costs function estimates the costs at every cycle.
Costs <- function (M_it) {
  # M_it: current health state
  # individual characteristic containing weight in third column
  c_it<-c()
  Merged_it <- cbind(df_x3,M_it,cost_i=0)
  colnames(Merged_it)<-c("ID","Sex","weight","sex_binary","M_it","cost_i")
  # update m_p_it with the appropriate probabilities
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="baseline" ,c_inf*weight+c_admin+c_monitor+c_mtxCo , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="ACR20" ,c_inf*weight+c_admin+c_monitor+c_mtxCo , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="ACR50" ,c_inf*weight+c_admin+c_monitor+c_mtxCo , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="ACR70" ,c_inf*weight+c_admin+c_monitor+c_mtxCo , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="nonresp" ,c_inf*weight+c_admin+c_monitor+c_mtxCo , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="Dead" ,0 , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="Discontinue" ,c_inf*weight+c_admin+c_monitor+c_ae_inf+c_mtxCo , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="last_line" ,c_mtx , cost_i))
  c_it<-Merged_it$cost_i
  return(c_it)  # return costs accrued this cycle
}

test_costs<-Costs(M_it=v_M_init)
#### 05.3 Health outcome function ####
# The Effs function to update the utilities at every cycle.

Effs <- function (M_it,HAQ_at,v_TS_at,v_last_dur_at,age_at) {
  # M_it: current health state
  q_it <- c()
  Merged_q_it<-cbind(df_x3,M_it,v_TS_at,HAQ_at,v_HAQ_init,v_last_dur_at,age_at,qaly=0)
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="baseline" , (.144)*(1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="ACR20" , (.144)*(1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="ACR50" , (.144)*(1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="ACR70" , (.144)*(1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="nonresp" , (.144)*(1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="discontinue" , (.144)*(1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))+u_ae_inf), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="last_line" & v_last_dur_at<15 , (.144)*((1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))+.0269*v_last_dur_at)), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="last_line" & v_last_dur_at>15 , (.144)*(1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))+.0269*15), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="Dead" , 0, qaly))
  q_it<-Merged_q_it$qaly
  return(q_it)  # return the QALYs accrued this cycle
}
test_Effs<-Effs(M_it=v_M_test_1,v_TS_at=v_Ts_init,HAQ_at=v_HAQ_init,v_last_dur_at=v_last_dur_init,age_at=v_age_init)

#### 06 Run Microsimulation ####
MicroSim <- function(n_i, seed = 1) {
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
  v_age_sim<-v_age_init#initialize age
  v_age_cont<-v_age_init#initialize cont_age
  v_last_dur<-v_last_dur_init
  v_last_dur_cont<-v_last_dur_init
  
  v_HAQ_at<-v_HAQ_init
  m_C[, 1] <- Costs(M_it=m_M[, 1])   # costs accrued individual i during cycle 0
  m_E[, 1] <- matrix(Effs(M_it=m_M[, 1],HAQ_at=v_HAQ_init,v_TS_at=v_Ts_init,v_last_dur_at=v_last_dur_init,age_at=v_age_sim))    # QALYs accrued individual i during cycle 0
  # open a loop for time running cycles 1 to n_t 
  for (t in 1:n_t) {
    v_p <- Probs(M_it=m_M[, t], v_haq=v_HAQ_at, v_age=v_age_sim, t)  # calculate the transition probabilities for the cycle based on health state t
    m_M[, t + 1]  <- samplev(v_p, 1)       # sample the current health state and store that state in matrix m_M 
    m_C[, t + 1]  <- Costs(M_it=m_M[, t + 1])   # calculate costs per individual during cycle t + 1
    m_E[, t + 1]  <- matrix(Effs (M_it=m_M[, t + 1],v_TS_at = v_Ts,HAQ_at =v_HAQ_at,v_last_dur_at = v_last_dur,age_at=v_age_sim))   # calculate QALYs per individual during cycle t + 1

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
  tc <- m_C%*%(v_dwc[1:(n_t+1)])     #total (discounted) cost per individual
  te <- m_E%*%(v_dwe[1:(n_t+1)])      #total (discounted) QALYs per individual 
  tc_hat <- mean(tc)     # average (discounted) cost 
  te_hat <- mean(te)     # average (discounted) QALYs
  
  # store the results from the simulation in a list
  results <- list(m_M = m_M, m_C = m_C, m_E = m_E, tc = tc , te = te, tc_hat = tc_hat, te_hat = te_hat)   
  return(results)  # return the results
} # end of the MicroSim function  

test_micro<-MicroSim(n_i = n_i)

#define functions for remicade

Probs2 <- function(M_it, v_haq, v_age_at,t) { 
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
  
  for (i in 1:n_i){
    m_p_it[1,i]<-0 #no one moves to baseline
    m_p_it[2,i]<-ifelse(M_it[i]=="baseline",(1-(1-p_mort[i])*p_discontinuer)*p_ACR20r,m_p_it[2,i])#moving to ACR20 from baseline
    m_p_it[3,i]<-ifelse(M_it[i]=="baseline",(1-(1-p_mort[i])*p_discontinuer)*p_ACR50r,m_p_it[3,i])#moving to ACR50 from baseline
    m_p_it[4,i]<-ifelse(M_it[i]=="baseline",(1-(1-p_mort[i])*p_discontinuer)*p_ACR70r,m_p_it[4,i])#moving to ACR50 from baseline
    m_p_it[5,i]<-ifelse(M_it[i]=="baseline",(1-(1-p_mort[i])*p_discontinuer)*p_nonrespr,m_p_it[5,i])#moving to nonresp from baseline
    m_p_it[6,i]<-ifelse(M_it[i]=="baseline",(1-p_mort[i])*p_discontinuer,m_p_it[6,i])#moving to discontinue from baseline
    m_p_it[7,i]<-ifelse(M_it[i]=="baseline",0,m_p_it[7,i])#moving to lastline from baseline
    m_p_it[8,i]<-ifelse(M_it[i]=="baseline",p_mort[i],m_p_it[8,i])#moving to death from baseline
    
    m_p_it[2,i]<-ifelse(M_it[i]=="ACR20",(1-(1-p_mort[i])*p_discontinuer)*p_ACR20r,m_p_it[2,i])
    m_p_it[3,i]<-ifelse(M_it[i]=="ACR20",(1-(1-p_mort[i])*p_discontinuer)*p_ACR50r,m_p_it[3,i])
    m_p_it[4,i]<-ifelse(M_it[i]=="ACR20",(1-(1-p_mort[i])*p_discontinuer)*p_ACR70r,m_p_it[4,i])
    m_p_it[5,i]<-ifelse(M_it[i]=="ACR20",(1-(1-p_mort[i])*p_discontinuer)*p_nonrespr,m_p_it[5,i])#moving to nonresp from baseline
    m_p_it[6,i]<-ifelse(M_it[i]=="ACR20",(1-p_mort[i])*p_discontinuer,m_p_it[6,i])#moving to discontinue from baseline
    m_p_it[7,i]<-ifelse(M_it[i]=="ACR20",0,m_p_it[7,i])#moving to lastline from baseline
    m_p_it[8,i]<-ifelse(M_it[i]=="ACR20",p_mort[i],m_p_it[8,i])#moving to death from baseline
    
    m_p_it[2,i]<-ifelse(M_it[i]=="ACR50",(1-(1-p_mort[i])*p_discontinuer)*p_ACR20r,m_p_it[2,i])
    m_p_it[3,i]<-ifelse(M_it[i]=="ACR50",(1-(1-p_mort[i])*p_discontinuer)*p_ACR50r,m_p_it[3,i])
    m_p_it[4,i]<-ifelse(M_it[i]=="ACR50",(1-(1-p_mort[i])*p_discontinuer)*p_ACR70r,m_p_it[4,i])
    m_p_it[5,i]<-ifelse(M_it[i]=="ACR50",(1-(1-p_mort[i])*p_discontinuer)*p_nonrespr,m_p_it[5,i])#moving to nonresp from baseline
    m_p_it[6,i]<-ifelse(M_it[i]=="ACR50",(1-p_mort[i])*p_discontinuer,m_p_it[6,i])#moving to discontinue from baseline
    m_p_it[7,i]<-ifelse(M_it[i]=="ACR50",0,m_p_it[7,i])#moving to lastline from baseline
    m_p_it[8,i]<-ifelse(M_it[i]=="ACR50",p_mort[i],m_p_it[8,i])#moving to death from baseline
    
    m_p_it[2,i]<-ifelse(M_it[i]=="ACR70",(1-(1-p_mort[i])*p_discontinuer)*p_ACR20r,m_p_it[2,i])
    m_p_it[3,i]<-ifelse(M_it[i]=="ACR70",(1-(1-p_mort[i])*p_discontinuer)*p_ACR50r,m_p_it[3,i])
    m_p_it[4,i]<-ifelse(M_it[i]=="ACR70",(1-(1-p_mort[i])*p_discontinuer)*p_ACR70r,m_p_it[4,i])
    m_p_it[5,i]<-ifelse(M_it[i]=="ACR70",(1-(1-p_mort[i])*p_discontinuer)*p_nonrespr,m_p_it[5,i])#moving to nonresp from baseline
    m_p_it[6,i]<-ifelse(M_it[i]=="ACR70",(1-p_mort[i])*p_discontinuer,m_p_it[6,i])#moving to discontinue from baseline
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

Costs2 <- function (M_it) {
  # M_it: current health state
  # individual characteristic containing weight in third column
  c_it<-c()
  Merged_it <- cbind(df_x3,M_it,cost_i=0)
  colnames(Merged_it)<-c("ID","Sex","weight","sex_binary","M_it","cost_i")
  # update m_p_it with the appropriate probabilities
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="baseline" ,c_rem*weight+c_admin+c_monitor+c_mtxCo , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="ACR20" ,c_rem*weight+c_admin+c_monitor+c_mtxCo , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="ACR50" ,c_rem*weight+c_admin+c_monitor+c_mtxCo , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="ACR70" ,c_rem*weight+c_admin+c_monitor+c_mtxCo , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="nonresp" ,c_rem*weight+c_admin+c_monitor+c_mtxCo , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="Dead" ,0 , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="Discontinue" ,c_rem*weight+c_admin+c_monitor+c_ae_rem+c_mtxCo , cost_i))
  Merged_it<-Merged_it %>% mutate(cost_i = ifelse(M_it=="last_line" ,c_mtx , cost_i))
  c_it<-Merged_it$cost_i
  return(c_it)  # return costs accrued this cycle
}

Effs2 <- function (M_it,HAQ_at,v_TS_at,v_last_dur_at,age_at) {
  # M_it: current health state
  q_it <- c()
  Merged_q_it<-cbind(df_x3,M_it,v_TS_at,HAQ_at,v_HAQ_init,v_last_dur_at,age_at,qaly=0)
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="baseline" , (.144)*(1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="ACR20" , (.144)*(1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="ACR50" , (.144)*(1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="ACR70" , (.144)*(1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="nonresp" , (.144)*(1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="discontinue" , (.144)*(1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at)))+u_ae_rem, qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="last_line" & v_last_dur_at<15 , (.144)*((1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))+.0269*v_last_dur_at)), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="last_line" & v_last_dur_at>15 , (.144)*((1-1/(1+exp(2.0734+0.0058*age_at+.0023*v_TS_at-.2004*v_HAQ_init-.2914*sex_binary+.0249*1-.8647*HAQ_at))+.0269*15)), qaly))
  Merged_q_it<-Merged_q_it %>% mutate(qaly = ifelse(M_it=="Dead" , 0, qaly))
  q_it<-Merged_q_it$qaly
  return(q_it)  # return the QALYs accrued this cycle
}

MicroSim2 <- function(n_i, seed = 1) {
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
  v_age_sim<-v_age_init #initialize age
  v_last_dur<-v_last_dur_init
  v_HAQ_at<-v_HAQ_init
  m_C[, 1] <- Costs2(M_it=m_M[, 1])   # costs accrued individual i during cycle 0
  m_E[, 1] <- matrix(Effs2(M_it=m_M[, 1],HAQ_at=v_HAQ_init,v_TS_at=v_Ts_init,v_last_dur_at=v_last_dur_init,age_at=v_age_sim))    # QALYs accrued individual i during cycle 0
  # open a loop for time running cycles 1 to n_t 
  for (t in 1:n_t) {
    v_p <- Probs2(M_it=m_M[, t], v_haq=v_HAQ_at, v_age=v_age_sim, t)  # calculate the transition probabilities for the cycle based on health state t
    m_M[, t + 1]  <- samplev(v_p, 1)       # sample the current health state and store that state in matrix m_M 
    m_C[, t + 1]  <- Costs2(M_it=m_M[, t + 1])   # calculate costs per individual during cycle t + 1
    m_E[, t + 1]  <- matrix(Effs2 (M_it=m_M[, t + 1],v_TS_at = v_Ts,HAQ_at =v_HAQ_at,v_last_dur_at = v_last_dur,age_at=v_age_sim))   # calculate QALYs per individual during cycle t + 1
    
    ifelse(v_age_sim<110,v_age_sim+1,v_age_sim)
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
    
    v_Ts<-if_else(m_M[, t + 1] != "Dead", v_Ts + 1, v_Ts) # update time since illness onset for t + 1 
    v_last_dur<-if_else(m_M[, t + 1] == "last_line", v_last_dur+1, v_last_dur)
    
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


# Run the simulation for both no treatment and treatment options
outcomes_inf  <- MicroSim(n_i)
outcomes_rem     <- MicroSim2(n_i)

#### 07 Visualize results ####
options(scipen = 999)

### infliximab
plot(density(outcomes_inf$tc), main = paste("INF Total cost per person"), xlab = "Cost ($)")
plot(density(outcomes_inf$te), main = paste("INF Total QALYs per person"), xlab = "QALYs")

n_t<-61
v_n   <- c("baseline","ACR20","ACR50","ACR70", "nonresp", "discontinue","last_line","Dead")          # vector with state names
n_s <- length(v_n)
plot_m_TR(outcomes_inf$m_M)    # health state trace
v_n   <- c("ACR20","ACR50","ACR70", "nonresp", "discontinue","last_line","Dead")          # vector with state names
n_t<-60
n_s <- length(v_n)

### remicade
plot(density(outcomes_rem$tc), main = paste("REM Total cost per person"), xlab = "Cost ($)")
plot(density(outcomes_rem$te), main = paste("REM Total QALYs per person"), xlab = "QALYs")

plot_m_TR(outcomes_rem$m_M)    # health state trace


#### 08 Cost Effectiveness Analysis ####
# store the mean costs of each strategy in a new variable C (vector of costs)
v_C <- c(outcomes_inf$tc_hat, outcomes_rem$tc_hat)
# store the mean QALYs of each strategy in a new variable E (vector of effects)
v_E <- c(outcomes_inf$te_hat, outcomes_rem$te_hat)

v_names_str <- c("INF", "REM") # strategy names

# use dampack to calculate the ICER
calculate_icers(cost       = v_C,
                effect     = v_E,
                strategies = v_names_str)





#################################
# 5.2 cost function : 2 ae's reported in the study --> Incidence of drug-related adverse events (35.2% vs 35.9%) and detection of antidrug antibodies (48.4% vs 48.2%) 
# cylce 7.5 weeks administered
#probabilities don't sum up to 1 
#Utility score relationship with HAQ, previous DMARD=0,1,2?
#where do hospital days, baseline missed work days, and unemployment fall in? cost?
################################