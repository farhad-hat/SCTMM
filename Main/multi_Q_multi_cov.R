# capture.output(sessionInfo(),file=file.path(pathScripts,"Rsession.txt" ))
library(tidyverse)
library(expm)
library(Metrics)
library(foreach)
library(doParallel)

transitions.table <- function(data){
  
  data$end.time <- NA
  data$end.time[-nrow(data)] <- data$time[2:nrow(data)]
  
  data$end.state <- data$state
  data$end.state[-nrow(data)] <- data$state[2:nrow(data)]
  
  data <- data %>% group_by(id) %>% slice(-n())
  
  data <- data %>% mutate(t = end.time - lag(end.time))
  data$t[is.na(data$t)] <- data$end.time[is.na(data$t)]
  
  data <- data.frame(data)
  
  return(data)
}


log_lik <- function(par, data){
  cat("Calculating likelihood ...", '\n')
  p <- ncol(data)-6
  Q.0 <- matrix(NA, nrow = nrow.Q, ncol = nrow.Q)
  r <- 1
  for (k in 1:nrow.Q) {
    for (l in c(1:nrow.Q)[c(1:nrow.Q)!=k]) {
      Q.0[k,l] <- par[r]
      r <- r+1
    }
    #Q.0[k,k] <- -sum(Q.0[k,-k])
  }
  
  
  beta <- list()
  for (s in 1:p) {
    beta[[s]] <- matrix(NA, nrow = nrow.Q, ncol = nrow.Q)
  }
  
  for (s in 1:p) {
    for (k in 1:nrow.Q) {
      for (l in c(1:nrow.Q)[c(1:nrow.Q)!=k]) {
        beta[[s]][k,l] <- par[r]
        r <- r+1
      }
    }
  }

  
  Q <- matrix(NA, nrow = nrow.Q, ncol = nrow.Q)
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  lLiki <- numeric(nrow(data))
  ll.foreach <- foreach(i=1:nrow(data), .combine="cbind",  .packages="expm",  .export = ls(globalenv())) %dopar% {
    #for(i in 1:nrow(data)){
    for (k in 1:nrow.Q) {
      for (l in c(1:nrow.Q)[c(1:nrow.Q)!=k]) {
        
        haz <- 0
        for (s in 1:p) {
          haz <- haz + (data[i,(s+2)]*beta[[s]][k,l])
        }
        
        Q[k,l] <- Q.0[k,l]*exp(haz)
      }
      Q[k,k] <- -sum(Q[k,-k])
    }
    
    Pts <- expm(Q*data$t[i])
    lLiki[i] <- log(Pts[data$state[i],
                        data$end.state[i]])
  }
  lLik <- sum(ll.foreach)
  #lLik <- sum(lLiki)
  stopCluster(cl)
  
  return(-lLik)  
  
  # return lLik (maximize)
  # return -2*lLik (minimize)
}


grad <- function(par, data){
  
  p <- ncol(data)-6
  Q.0 <- matrix(NA, nrow = nrow.Q, ncol = nrow.Q)
  r <- 1
  for (k in 1:nrow.Q) {
    for (l in c(1:nrow.Q)[c(1:nrow.Q)!=k]) {
      Q.0[k,l] <- par[r]
      r <- r+1
    }
    Q.0[k,k] <- -sum(Q.0[k,-k])
  }
  
  
  beta <- list()
  for (s in 1:p) {
    beta[[s]] <- matrix(NA, nrow = nrow.Q, ncol = nrow.Q)
  }
  
  for (s in 1:p) {
    for (k in 1:nrow.Q) {
      for (l in c(1:nrow.Q)[c(1:nrow.Q)!=k]) {
        beta[[s]][k,l] <- par[r]
        r <- r+1
      }
    }
  }
  
  gr_sum <- numeric((nrow.Q*(nrow.Q-1))+(p*nrow.Q*(nrow.Q-1)))
  ####### wrt Q.0
  #r <- 1
  
  cat("Calculating derivatives wrt Q.0 ...", '\n')
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  gr_sum_foreach <- matrix(NA,nrow.Q,nrow.Q)
  gr.foreach <- foreach(n=1:nrow.Q, .combine='cbind') %:%
    foreach(m=c(1:nrow.Q)[c(1:nrow.Q)!=n], .combine='c', .packages='expm', .export = ls(globalenv())) %dopar% {
      
      #for (n in 1:nrow.Q) {
      
      #for (m in c(1:nrow.Q)[c(1:nrow.Q)!=n]) {
      
      der <- list()
      for(i in 1:nrow(data)){
        
        Q <- matrix(NA, nrow = nrow.Q, ncol = nrow.Q)
        for (k in 1:nrow.Q) {
          for (l in c(1:nrow.Q)[c(1:nrow.Q)!=k]) {
            
            haz <- 0
            for (s in 1:p) {
              haz <- haz + (data[i,(s+2)]*beta[[s]][k,l])
            }
            
            Q[k,l] <- Q.0[k,l]*exp(haz)
          }
          Q[k,k] <- -sum(Q[k,-k])
        }
        
        B <- Q*data$t[i]
        C <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
        
        haz <- 0
        for (s in 1:p) {
          haz <- haz + (data[i,(s+2)]*beta[[s]][n,m])
        }
        
        C[n,m] <- exp(haz) * data$t[i]
        H <- rbind(Q*data$t[i],matrix(0,nrow.Q,nrow.Q))
        H <- cbind(H, rbind(C,Q*data$t[i]))
        H <- expm(H)[1:nrow.Q,((nrow.Q+1):(2*nrow.Q))]
        der[[i]] <- H / expm(Q*data$t[i])
      }
      # Likelihood
      gr <- numeric(nrow(data))
      for(i in 1:nrow(data)){
        gr[i] <- der[[i]][data$state[i],
                          data$end.state[i]]
      }
      #gr_sum[r] <- sum(gr)
      #average relaxes any mini-batch size
      gr_sum_foreach[n,m] <- sum(gr) / nrow(data)
      #r <- r+1
    }
  gr_sum[1:(nrow.Q*(nrow.Q-1))] <- c(gr.foreach)
  #print(n)
  #}
  stopCluster(cl)
  ####### wrt beta
  cat("Calculating derivatives wrt beta ...", '\n')
  for (u in 1:p) {
    
    cores <- detectCores()
    cl <- makeCluster(cores[1]-1) #not to overload your computer
    registerDoParallel(cl)
    
    gr_sum_foreach <- matrix(NA,nrow.Q,nrow.Q)
    gr.foreach <- foreach(n=1:nrow.Q, .combine='cbind') %:%
      foreach(m=c(1:nrow.Q)[c(1:nrow.Q)!=n], .combine='c', .packages='expm', .export = ls(globalenv())) %dopar% {
        
        #for (n in 1:nrow.Q) {
        #for (m in c(1:nrow.Q)[c(1:nrow.Q)!=n]) {
        
        der <- list()
        for(i in 1:nrow(data)){
          
          Q <- matrix(NA, nrow = nrow.Q, ncol = nrow.Q)
          for (k in 1:nrow.Q) {
            for (l in c(1:nrow.Q)[c(1:nrow.Q)!=k]) {
              
              haz <- 0
              for (s in 1:p) {
                haz <- haz + (data[i,(s+2)]*beta[[s]][k,l])
              }
              
              Q[k,l] <- Q.0[k,l]*exp(haz)
            }
            Q[k,k] <- -sum(Q[k,-k])
          }
          
          B <- Q*data$t[i]
          C <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
          
          haz <- 0
          for (s in 1:p) {
            haz <- haz + (data[i,(s+2)]*beta[[s]][n,m])
          }
          
          C[n,m] <- (Q.0[n,m]*data[i,(u+2)]*exp(haz)) * data$t[i]
          H <- rbind(Q*data$t[i],matrix(0,nrow.Q,nrow.Q))
          H <- cbind(H, rbind(C,Q*data$t[i]))
          H <- expm(H)[1:nrow.Q,((nrow.Q+1):(2*nrow.Q))]
          der[[i]] <- H / expm(Q*data$t[i])
        }
        # Likelihood
        gr <- numeric(nrow(data))
        for(i in 1:nrow(data)){
          gr[i] <- der[[i]][data$state[i],
                            data$end.state[i]]
        }
        #gr_sum[r] <- sum(gr)
        gr_sum_foreach[n,m] <- sum(gr) / nrow(data)
        #r <- r+1
      }
    gr_sum[ ((u*nrow.Q*(nrow.Q-1))+1):((u+1)*nrow.Q*(nrow.Q-1)) ] <- c(gr.foreach)
    stopCluster(cl)
    #}
  }
  return(-gr_sum)  
}


stc_grad_desc <- function(init.par, data, batch.size, max.iter, lambda, tol){
  ptm <- proc.time()
  batch.ids.from.all.states <- vector()
  for (state in sort(unique(data$state))) {
    data.temp <- data[(data$state == state) & (!(data$id %in% batch.ids.from.all.states)),]
    batch.ids <- sample(unique(data.temp$id), 1)
    batch.ids.from.all.states <- append(batch.ids.from.all.states, batch.ids)
  }
  data.temp <- data[!(data$id %in% batch.ids.from.all.states),]
  batch.ids <- sample(unique(data.temp$id), batch.size)
  batch.ids.from.all.states <- append(batch.ids.from.all.states, batch.ids)
  batch.data <- data[data$id %in% batch.ids.from.all.states,]
  
  par <- list()
  par[[1]] <- init.par
  ll <- vector()
  init.lambda <- lambda
  
  gradient.old <- grad(par[[1]], batch.data)
  ll[1] <- log_lik(par=par[[1]], data=data)
  ll[2] <- log_lik(par=(par[[1]] - lambda*gradient.old), data=data)
  
  stop = FALSE
  for(step in 2:max.iter) {
    lambda <- init.lambda
    
    batch.ids.from.all.states <- vector()
    for (state in sort(unique(data$state))) {
      data.temp <- data[(data$state == state) & (!(data$id %in% batch.ids.from.all.states)),]
      batch.ids <- sample(unique(data.temp$id), 1)
      batch.ids.from.all.states <- append(batch.ids.from.all.states, batch.ids)
    }
    data.temp <- data[!(data$id %in% batch.ids.from.all.states),]
    batch.ids <- sample(unique(data.temp$id), batch.size)
    batch.ids.from.all.states <- append(batch.ids.from.all.states, batch.ids)
    batch.data <- data[data$id %in% batch.ids.from.all.states,]
    
    if(step>2){
      gradient.old <- grad(par[[step-1]], batch.data)
      #ll[step-1] <- log_lik(par=par[[step-1]], data=data)
      ll[step] <- log_lik(par=(par[[step-1]] - lambda*gradient.old), data=data)
    }
    cat("Likelihood =",ll[step-1], '\n')
    while(ll[step] > ll[step-1]){
      lambda <- lambda/1.6
      if((lambda < 1e-10) || (ll[step] == -Inf)) {
        stop = TRUE
        break
      }
      cat("New learing rate =",lambda, '\n')
      cat("Likelihood =",ll[step], '\n')
      if(stop==FALSE) ll[step] <- log_lik(par=(par[[step-1]] - lambda*gradient.old), data=data)
    }
    
    par[[step]] <- par[[step-1]] - lambda*gradient.old
    
    if(stop) break
    
    if(abs(ll[step]-ll[step-1]) < tol) break
    cat("ll = ", ll[step], '\n')
    cat("step = ", step, '\n')
  }
  ptm <- proc.time() - ptm
  return(list(par, ll, ptm))
}


# msm <- msm(formula = state ~ time,
#            subject = id,
#            data = data,
#            #qmatrix = rbind(c(-0.2, 0.1, 0.1),c(0.1, -0.2, 0.1),c(0.1, 0.1, -0.2)),
#            qmatrix = matrix(c(-3.5,rep(0.5,7),  
#                               rep(0.5,1),-3.5,rep(0.5,6),
#                               rep(0.5,2),-3.5,rep(0.5,5),
#                               rep(0.5,3),-3.5,rep(0.5,4),
#                               rep(0.5,4),-3.5,rep(0.5,3),
#                               rep(0.5,5),-3.5,rep(0.5,2),
#                               rep(0.5,6),-3.5,rep(0.5,1),
#                               rep(0.5,7),-3.5), byrow = T, nrow = nrow.Q),
#            covariates = ~cov_1+cov_2,
#            control = list(reltol = 1e-8,fnscale=4000, maxit = 10000,trace=1,REPORT=1),
#            exacttimes = TRUE)
# 
# 
# qmatrix.msm(msm, ci='none')
# haz <- hazard.msm(msm)
# haz$cov_1




# 
# 
# library(optimx)
# hess <- gHgen(par = c(1.93884192, 0.61329431, 1.34655721, 3.55947315, 1.17057902, 3.01718844,
#                       0.14186130, 0.04599288, 0.09866043, 0.27111451, 0.09107962, 0.24998815,
#                       0.15324092, 0.04344108, 0.10835572, 0.28629542, 0.09421699, 0.22412849),
#               fn=log_lik, gr=grad, data=data)
# 
# ##########
# optim.res2 <- optim(par = c(0.1,0.1,0.1,0.1,0.1,0.1, 
#                             0.001,0.001,0.001,0.001,0.001,0.001,
#                             #0.001,0.001,0.001,0.001,0.001,0.001,
#                             #0.001,0.001,0.001,0.001,0.001,0.001,
#                             0.001,0.001,0.001,0.001,0.001,0.001), 
#                     fn = log_lik, gr=grad, data=data, method = 'BFGS', hessian = T,
#                     control=list(reltol = 1e-5,fnscale=4000,maxit = 1000,trace=1,REPORT=1))
# 
# 
# fisher_info<-solve(res$hessian)
#if -2*log.lik
#fisher_info<-solve(0.5*res$hessian)
# prop_sigma<-sqrt(diag(fisher_info))
# prop_sigma<-diag(prop_sigma)
# upper<-res$par+1.959964*prop_sigma
# lower<-res$par-1.959964*prop_sigma
# interval<-data.frame(value=res$par, upper=upper, lower=lower)
# 
# cbind(diag(lower),
#       diag(upper))
# 
# saveRDS(optim.res2, 'optim.res2.RData')
# 
# # msm <- msm(formula = state ~ time,
# #            subject = id,
# #            data = data,
# #            qmatrix = rbind(c(0, 0.4),c(0.4, 0)),
# #            covariates = ~cov_1,
# #            control = list(reltol = 1e-8,fnscale=4000, maxit = 1000,trace=1,REPORT=1))
# # 
# # 
# # qmatrix.msm(msm, ci='none')
# # pmatrix.msm(msm, t=9)
# 
# 
# 
# 
# data <- data.frame(CMCSimulation(Q,0,10))
# colnames(data) <- c('state','time')
# data$state <- data$state + 1
# data$id <- c(1)
# for (id in 2:4) {
#   temp <- data.frame(CMCSimulation(Q,0,10))
#   colnames(temp) <- c('state','time')
#   temp$state <- temp$state + 1
#   temp$id <- c(id)
#   
#   data <- rbind(data, temp)
# }
# 
# 
