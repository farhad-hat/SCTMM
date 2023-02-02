library(castor) 
library(tidyverse)
library(zoo)

 # Q.0 <- matrix(c(-0.5,0.4,0.1,  0.25,-1,0.75,  0.25,0.75,-1), byrow = T, nrow = 3)
 # beta <- c(0.1,0.1)
 # p <- 2

# Simulate 1 transition ---------------------------------------------------

initial.state <- function(Q.0, beta, p){
  z <- vector()
  for (r in 1:p) {
    z[r] <- 1 
  }
  
  
  Q <- matrix(NA, nrow = nrow(Q.0), ncol = nrow(Q.0))
  for (k in 1:nrow(Q)) {
    for (l in c(1:nrow(Q))[c(1:nrow(Q))!=k]) {
      
      haz <- 0
      for (r in 1:p) {
        haz <- haz + (z[r]*beta[r])
      }
      
      Q[k,l] <- Q.0[k,l]*exp(haz)
    }
    Q[k,k] <- -sum(Q[k,-k])
  }
  
  stationary <- sample(x = 1:(nrow(Q)), size = 1, prob = get_stationary_distribution(Q))
  return(list('stationary'=stationary,'z'=z))
} 

#initial.state(Q.0, beta, p)
#########################################################################################
simulate.transition <- function(Q.0, beta, current.state, z, p){
  
  Q <- matrix(NA, nrow = nrow(Q.0), ncol = nrow(Q.0))
  for (k in 1:nrow(Q)) {
    for (l in c(1:nrow(Q))[c(1:nrow(Q))!=k]) {
      
      haz <- 0
      for (r in 1:p) {
        haz <- haz + (z[r]*beta[r])
      }
      
      Q[k,l] <- Q.0[k,l]*exp(haz)
    }
    Q[k,k] <- -sum(Q[k,-k])
  }
  
  t = rexp(n=1, rate = -(Q[current.state, current.state]))
  s = sample(x = (1:nrow(Q))[-current.state], size = 1,
             prob = Q[current.state, -current.state] / -Q[current.state, current.state]) 
  
  return(c(s, t, z))
}

#simulate.transition(Q.0, beta, current.state=initial.state(Q.0, beta, p)$stationary, z=initial.state(Q.0, beta, p)$z, p)
#########################################################################################

simulate.id <- function(Q.0, beta, followup, id, p){
  
  surv.table <- matrix(c(initial.state(Q.0, beta, p)$stationary, 
                         0,
                         initial.state(Q.0, beta, p)$z), byrow = F, ncol = p+2)
  colnames(surv.table) <- c('state', 'time', paste("cov_", 1:p, sep = ""))
  surv.table <- data.frame(surv.table)
  
  while(surv.table$time[nrow(surv.table)] < followup){
    new.obs <- simulate.transition(Q.0=Q.0, beta=beta, current.state=surv.table[nrow(surv.table), 1], z=abs(rnorm(p)), p) + c(0, surv.table[nrow(surv.table), 2], rep(0,p))
    
    if(new.obs[2] > followup) break
    
    floor <- floor((new.obs[2] - surv.table[nrow(surv.table), 2]) / 0.25)
    if(floor != 0){
      for(i in 1:floor){
        surv.table <- rbind(surv.table, c(surv.table[,1][dim(surv.table)[1]], 
                                          (floor(surv.table[,2][dim(surv.table)[1]]*4)/4) + (0.25),
                                          unname(surv.table[,3:dim(surv.table)[2]][dim(surv.table)[1],])))
      }
    }
    surv.table <- rbind(surv.table, new.obs)
  }
  
  surv.table$id <- rep(id,nrow(surv.table))
  
  return(surv.table)
}

#simulate.id(Q.0=Q.0, beta=beta, followup=3, id=1, p)
##########################################################################################


# Simulate Multiple Observations ------------------------------------------

# N=10;current.state = 1;followup=5;Q = matrix(c(-0.5,0.4,0.1,  0.25,-1,0.75,  0,0,0),byrow = T,nrow=3)
simulate.data <- function(Q.0, beta, followup, N, p){
  
  sim.data <- simulate.id(Q.0=Q.0, beta=beta, followup=followup, id=1, p=p)
  
  for(id in 2:N){
    sim.data <- rbind(sim.data, 
                      simulate.id(Q.0=Q.0, beta=beta, followup=followup, id=id, p=p))  
  }
  
  #sim.data$state <- factor(sim.data$state)
  # sim.data <- sim.data %>% group_by(id) %>% mutate(T = lead(time) - time)
  # sim.data$T[is.na(sim.data$T)] <- followup - sim.data$time[is.na(sim.data$T)]
  # sim.data <- as.data.frame(sim.data)
  # 
  # sim.data <- sim.data %>% group_by(id) %>% mutate(end.state = lead(state))
  # sim.data <- data.frame(sim.data)
  # sim.data$end.state <- na.locf(sim.data$end.state)
  return(sim.data)  
}

#simulate.data(Q.0=Q.0, beta=beta, followup=15, N=10, p)
