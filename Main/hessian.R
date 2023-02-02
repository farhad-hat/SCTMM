library(foreach)
library(doParallel)
#library(MASS)

hess.function <- function(par, data){
  
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
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  hes.list <- list()
  sum.hessian <- 0
  #for(mk in 1:nrow(data)){
  ll.foreach <- foreach(mk=1:nrow(data)) %dopar% {
    # Calculating hessian
    # Calculating second derivatives of Q
    
    Q <- matrix(NA, nrow = nrow.Q, ncol = nrow.Q)
    for (k in 1:nrow.Q) {
      for (l in c(1:nrow.Q)[c(1:nrow.Q)!=k]) {
        
        haz <- 0
        for (s in 1:p) {
          haz <- haz + (data[mk,(s+2)]*beta[[s]][k,l])
        }
        
        Q[k,l] <- Q.0[k,l]*exp(haz)
      }
    }
    
    
    state <- data$state[mk]
    end.state <- data$end.state[mk]
    
    haz <- 0
    for (s in 1:p) {
      haz <- haz + (data[mk,(s+2)]*beta[[s]][state, end.state])
    }
    haz <- exp(haz)
    
    # d.q.Q <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
    # d.q.Q[state, end.state] <- haz
    # diag(d.q.Q) <- NA
    
    
    # d.beta.Q <- list()
    # for (s in 1:p) {
    #   d.beta.Q[[s]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
    #   d.beta.Q[[s]][state, end.state] <- Q.0[state, end.state] * data[mk,(s+2)] * haz
    #   diag(d.beta.Q[[s]]) <- NA
    # }
    
    
    d.q.q.Q <- list()
    d.beta.q.Q <- list(list())
    for (s in 1:(nrow.Q*(nrow.Q-1))) {
      d.q.q.Q[[s]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
      diag(d.q.q.Q[[s]]) <- NA
      
      
      d.beta.q.Q[[s]] <- list()
      for (r in 1:p) {
        d.beta.q.Q[[s]][[r]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
        d.beta.q.Q[[s]][[r]][state, end.state] <- data[mk,(r+2)] * haz
        diag(d.beta.q.Q[[s]][[r]]) <- NA
      }
    }
    
    
    d.q.beta.Q <- list(list())
    d.beta.beta.Q <- list(list(list()))
    for (s in 1:p) {
      d.q.beta.Q[[s]] <- list()
      
      
      d.beta.beta.Q[[s]] <- list(list())
      
      
      
      for (r in 1:(nrow.Q*(nrow.Q-1))) {
        d.q.beta.Q[[s]][[r]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
        d.q.beta.Q[[s]][[r]][state, end.state] <- data[mk,(s+2)] * haz
        diag(d.q.beta.Q[[s]][[r]]) <- NA
        
        
        d.beta.beta.Q[[s]][[r]] <- list()
        for (u in 1:p) {
          d.beta.beta.Q[[s]][[r]][[u]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
          d.beta.beta.Q[[s]][[r]][[u]][state, end.state] <-  data[mk,(s+2)]^2 * Q.0[state, end.state] * haz
          diag(d.beta.beta.Q[[s]][[r]][[u]]) <- NA
        }
      }
    }
    
    
    
    # vec.d.q.Q <- as.vector(t(d.q.Q))  
    # vec.d.q.Q <- vec.d.q.Q[!is.na(vec.d.q.Q)]
    # 
    # 
    # vec.d.beta.Q <- as.vector(t(do.call(rbind, d.beta.Q)))
    # vec.d.beta.Q <- vec.d.beta.Q[!is.na(vec.d.beta.Q)]
    
    
    vec.d.q.q.Q <- as.vector(t(do.call(rbind, d.q.q.Q)))
    vec.d.q.q.Q <- vec.d.q.q.Q[!is.na(vec.d.q.q.Q)]
    
    
    vec.d.beta.q.Q <- vector()
    for (s in 1:(nrow.Q*(nrow.Q-1))) {
      vec.d.beta.q.Q <- append(vec.d.beta.q.Q, t(do.call(rbind, d.beta.q.Q[[s]])))
    }
    vec.d.beta.q.Q <- vec.d.beta.q.Q[!is.na(vec.d.beta.q.Q)]
    
    
    vec.d.q.beta.Q <- vector()
    vec.d.beta.beta.Q <- vector()
    for (s in 1:p) {
      vec.d.q.beta.Q <- append(vec.d.q.beta.Q, t(do.call(rbind, d.q.beta.Q[[s]])))
      
      
      for (r in 1:(nrow.Q*(nrow.Q-1))) {
        vec.d.beta.beta.Q <- append(vec.d.beta.beta.Q, t(do.call(rbind, d.beta.beta.Q[[s]][[r]])))
      }
    }
    vec.d.q.beta.Q <- vec.d.q.beta.Q[!is.na(vec.d.q.beta.Q)]
    
    
    
    vec.d.beta.beta.Q <- vec.d.beta.beta.Q[!is.na(vec.d.beta.beta.Q)]
    
    
    hessian.up <- cbind(matrix(vec.d.q.q.Q, nrow=nrow.Q*(nrow.Q-1), ncol=nrow.Q*(nrow.Q-1), byrow=TRUE), 
                        matrix(vec.d.beta.q.Q, nrow=nrow.Q*(nrow.Q-1), ncol=p*(nrow.Q*(nrow.Q-1)), byrow=TRUE))
    hessian.down <- cbind(matrix(vec.d.q.beta.Q, nrow=p*(nrow.Q*(nrow.Q-1)), ncol=nrow.Q*(nrow.Q-1), byrow=TRUE), 
                          matrix(vec.d.beta.beta.Q, nrow=p*(nrow.Q*(nrow.Q-1)), ncol=p*(nrow.Q*(nrow.Q-1)), byrow=TRUE))
    hessian.part.Q <- rbind(hessian.up, hessian.down)  
    
    # Calculating second derivatives of Q^2
    ####################################### First situation
    # 1
    d.q.q.Q2 <- list()
    for (s in 1:(nrow.Q*(nrow.Q-1))) {
      d.q.q.Q2[[s]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
      diag(d.q.q.Q2[[s]]) <- NA
    }
    
    s <- 1
    ind.d.q.q.Q2 <- data.frame()
    for (a in 1:nrow.Q) {
      for (b in c(1:nrow.Q)[c(1:nrow.Q)!=a]) {
        ind.d.q.q.Q2 <- rbind(ind.d.q.q.Q2, c(s,a,b))
        s <- s+1
      }
    }
    colnames(ind.d.q.q.Q2) <- c('s','a','b')
    
    s <- 1
    for (i in 1:nrow.Q) {
      for (j in c(1:nrow.Q)[c(1:nrow.Q)!=i]) {
        a <- ind.d.q.q.Q2[ind.d.q.q.Q2$s==s,]$a
        b <- ind.d.q.q.Q2[ind.d.q.q.Q2$s==s,]$b
        
        if(i==state && j!=end.state) {
          
          if(a==j && b==end.state){
            haz.ab <- 0
            for (u in 1:p) {
              haz.ab <- haz.ab + (data[mk,(u+2)]*beta[[u]][a, b])
            }
            haz.ab <- exp(haz.ab)
            
            haz.state.j <- 0
            for (u in 1:p) {
              haz.state.j <- haz.state.j + (data[mk,(u+2)]*beta[[u]][state, j])
            }
            haz.state.j <- exp(haz.state.j)
            
            d.q.q.Q2[[s]][a,b] <- haz.ab * haz.state.j
          }
        }
        s <- s+1
      }
    }
    
    # 2
    d.beta.q.Q2 <- list(list())
    for (s in 1:(nrow.Q*(nrow.Q-1))) {
      d.beta.q.Q2[[s]] <- list()
    }
    for (s in 1:(nrow.Q*(nrow.Q-1))) {
      for (r in 1:p) {
        d.beta.q.Q2[[s]][[r]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
        diag(d.beta.q.Q2[[s]][[r]]) <- NA
      }
    }
    
    s <- 1
    ind.d.beta.q.Q2 <- data.frame()
    for (c in 1:p) {
      for (a in 1:nrow.Q) {
        for (b in c(1:nrow.Q)[c(1:nrow.Q)!=a]) {
          ind.d.beta.q.Q2 <- rbind(ind.d.beta.q.Q2, c(s,a,b,c))
          s <- s+1
        }
      }
    }
    colnames(ind.d.beta.q.Q2) <- c('s','a','b','c')
    
    
    s <- 1
    for (i in 1:nrow.Q) {
      for (j in c(1:nrow.Q)[c(1:nrow.Q)!=i]) {
        a <- ind.d.beta.q.Q2[ind.d.beta.q.Q2$s==s,]$a
        b <- ind.d.beta.q.Q2[ind.d.beta.q.Q2$s==s,]$b
        c <- ind.d.beta.q.Q2[ind.d.beta.q.Q2$s==s,]$c
        
        if(i==state && j!=end.state) {
          
          if((a==j && b==end.state) ||  (a==state && b==j)) {
            haz.state.j <- 0
            for (u in 1:p) {
              haz.state.j <- haz.state.j + (data[mk,(u+2)]*beta[[u]][state, j])
            }
            haz.state.j <- exp(haz.state.j)
            
            d.beta.q.Q2[[s]][[c]][a,b] <- Q[j,end.state] * data[mk,(c+2)] * haz.state.j
          }
        }
        s <- s+1
      }
    }
    
    # 3
    d.q.beta.Q2 <- list(list())
    for (s in 1:p) {
      d.q.beta.Q2[[s]] <- list()
    }
    for (s in 1:p) {
      for (r in 1:(nrow.Q*(nrow.Q-1))) {
        d.q.beta.Q2[[s]][[r]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
        diag(d.q.beta.Q2[[s]][[r]]) <- NA
      }
    }
    
    s <- 1
    ind.d.q.beta.Q2 <- data.frame()
    for (a in 1:nrow.Q) {
      for (b in c(1:nrow.Q)[c(1:nrow.Q)!=a]) {
        for (c in 1:p) {
          ind.d.q.beta.Q2 <- rbind(ind.d.q.beta.Q2, c(s,a,b,c))
          s <- s+1
        }
      }
    }
    colnames(ind.d.q.beta.Q2) <- c('s','a','b','c')
    
    
    s <- 1
    for (i in 1:nrow.Q) {
      for (j in c(1:nrow.Q)[c(1:nrow.Q)!=i]) {
        a <- ind.d.q.beta.Q2[ind.d.q.beta.Q2$s==s,]$a
        b <- ind.d.q.beta.Q2[ind.d.q.beta.Q2$s==s,]$b
        c <- ind.d.q.beta.Q2[ind.d.q.beta.Q2$s==s,]$c
        
        if(i==state && j!=end.state) {
          
          if((a==state && b==j)) {
            haz.ab <- 0
            for (u in 1:p) {
              haz.ab <- haz.ab + (data[mk,(u+2)]*beta[[u]][a, b])
            }
            haz.ab <- exp(haz.ab)
            
            d.q.beta.Q2[[c]][[s]][a,b] <- haz.ab * data[mk,(c+2)] * Q[j,end.state] 
          }
          
          if((a==j && b==end.state)) {
            haz.ab <- 0
            for (u in 1:p) {
              haz.ab <- haz.ab + (data[mk,(u+2)]*beta[[u]][a, b])
            }
            haz.ab <- exp(haz.ab)
            d.q.beta.Q2[[c]][[s]][a,b] <- Q[state,j] * data[mk,(c+2)] * haz.ab 
          }
        }
        s <- s+1
      }
    }
    
    # 4
    d.beta.beta.Q2 <- list(list(list()))
    for (s in 1:p) {
      d.beta.beta.Q2[[s]] <- list(list())
    }
    for (s in 1:p) {
      for (r in 1:(nrow.Q*(nrow.Q-1))) {
        d.beta.beta.Q2[[s]][[r]] <- list()
      }
    }
    for (s in 1:p) {
      for (r in 1:(nrow.Q*(nrow.Q-1))) {
        for (u in 1:p) {
          d.beta.beta.Q2[[s]][[r]][[u]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
          diag(d.beta.beta.Q2[[s]][[r]][[u]]) <- NA
        }
      }
    }
    
    s <- 1
    ind.d.beta.beta.Q2 <- data.frame()
    for (c in 1:p) {
      for (r in 1:p) {
        for (a in 1:nrow.Q) {
          for (b in c(1:nrow.Q)[c(1:nrow.Q)!=a]) {
            
            ind.d.beta.beta.Q2 <- rbind(ind.d.beta.beta.Q2, c(s,a,b,c,r))
            s <- s+1
          }
        }
      }
    }
    colnames(ind.d.beta.beta.Q2) <- c('s','a','b','c','r')
    
    
    s <- 1
    for (i in 1:nrow.Q) {
      for (j in c(1:nrow.Q)[c(1:nrow.Q)!=i]) {
        a <- ind.d.beta.beta.Q2[ind.d.beta.beta.Q2$s==s,]$a
        b <- ind.d.beta.beta.Q2[ind.d.beta.beta.Q2$s==s,]$b
        c <- ind.d.beta.beta.Q2[ind.d.beta.beta.Q2$s==s,]$c
        r <- ind.d.beta.beta.Q2[ind.d.beta.beta.Q2$s==s,]$r
        
        if(i==state && j!=end.state) {
          
          if((a==state && b==j) ||  (a==j && b==end.state)) {
            d.beta.beta.Q2[[c]][[s]][[r]][a,b] <- Q[state,j] * data[mk,(r+2)] * data[mk,(c+2)] * Q[j,end.state]
          }
          
        }
        s <- s+1
      }
    }
    ####################################### 
    ####################################### Second situation
    # 1
    d.q.q.Q2 <- list()
    for (s in 1:(nrow.Q*(nrow.Q-1))) {
      d.q.q.Q2[[s]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
      diag(d.q.q.Q2[[s]]) <- NA
    }
    
    s <- 1
    ind.d.q.q.Q2 <- data.frame()
    for (a in 1:nrow.Q) {
      for (b in c(1:nrow.Q)[c(1:nrow.Q)!=a]) {
        ind.d.q.q.Q2 <- rbind(ind.d.q.q.Q2, c(s,a,b))
        s <- s+1
      }
    }
    colnames(ind.d.q.q.Q2) <- c('s','a','b')
    
    s <- 1
    for (i in 1:nrow.Q) {
      for (j in c(1:nrow.Q)[c(1:nrow.Q)!=i]) {
        a <- ind.d.q.q.Q2[ind.d.q.q.Q2$s==s,]$a
        b <- ind.d.q.q.Q2[ind.d.q.q.Q2$s==s,]$b
        
        if(i!=state && j==end.state) {
          
          if(a==state && b==i){
            haz.ab <- 0
            for (u in 1:p) {
              haz.ab <- haz.ab + (data[mk,(u+2)]*beta[[u]][a, b])
            }
            haz.ab <- exp(haz.ab)
            
            haz.i.end.state <- 0
            for (u in 1:p) {
              haz.i.end.state <- haz.i.end.state + (data[mk,(u+2)]*beta[[u]][i, end.state])
            }
            haz.i.end.state <- exp(haz.i.end.state)
            
            d.q.q.Q2[[s]][a,b] <- haz.ab * haz.state.j
          }
        }
        s <- s+1
      }
    }
    
    # 2
    d.beta.q.Q2 <- list(list())
    for (s in 1:(nrow.Q*(nrow.Q-1))) {
      d.beta.q.Q2[[s]] <- list()
    }
    for (s in 1:(nrow.Q*(nrow.Q-1))) {
      for (r in 1:p) {
        d.beta.q.Q2[[s]][[r]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
        diag(d.beta.q.Q2[[s]][[r]]) <- NA
      }
    }
    
    s <- 1
    ind.d.beta.q.Q2 <- data.frame()
    for (c in 1:p) {
      for (a in 1:nrow.Q) {
        for (b in c(1:nrow.Q)[c(1:nrow.Q)!=a]) {
          ind.d.beta.q.Q2 <- rbind(ind.d.beta.q.Q2, c(s,a,b,c))
          s <- s+1
        }
      }
    }
    colnames(ind.d.beta.q.Q2) <- c('s','a','b','c')
    
    
    s <- 1
    for (i in 1:nrow.Q) {
      for (j in c(1:nrow.Q)[c(1:nrow.Q)!=i]) {
        a <- ind.d.beta.q.Q2[ind.d.beta.q.Q2$s==s,]$a
        b <- ind.d.beta.q.Q2[ind.d.beta.q.Q2$s==s,]$b
        c <- ind.d.beta.q.Q2[ind.d.beta.q.Q2$s==s,]$c
        
        if(i!=state && j==end.state) {
          
          if((a==state && b==i) ||  (a==i && b==end.state)) {
            haz.i.end.state <- 0
            for (u in 1:p) {
              haz.i.end.state <- haz.i.end.state + (data[mk,(u+2)]*beta[[u]][i, end.state])
            }
            haz.i.end.state <- exp(haz.i.end.state)
            
            d.beta.q.Q2[[s]][[c]][a,b] <- Q[state,i] * data[mk,(c+2)] * haz.i.end.state
          }
        }
        s <- s+1
      }
    }
    
    # 3
    d.q.beta.Q2 <- list(list())
    for (s in 1:p) {
      d.q.beta.Q2[[s]] <- list()
    }
    for (s in 1:p) {
      for (r in 1:(nrow.Q*(nrow.Q-1))) {
        d.q.beta.Q2[[s]][[r]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
        diag(d.q.beta.Q2[[s]][[r]]) <- NA
      }
    }
    
    s <- 1
    ind.d.q.beta.Q2 <- data.frame()
    for (a in 1:nrow.Q) {
      for (b in c(1:nrow.Q)[c(1:nrow.Q)!=a]) {
        for (c in 1:p) {
          ind.d.q.beta.Q2 <- rbind(ind.d.q.beta.Q2, c(s,a,b,c))
          s <- s+1
        }
      }
    }
    colnames(ind.d.q.beta.Q2) <- c('s','a','b','c')
    
    
    s <- 1
    for (i in 1:nrow.Q) {
      for (j in c(1:nrow.Q)[c(1:nrow.Q)!=i]) {
        a <- ind.d.q.beta.Q2[ind.d.q.beta.Q2$s==s,]$a
        b <- ind.d.q.beta.Q2[ind.d.q.beta.Q2$s==s,]$b
        c <- ind.d.q.beta.Q2[ind.d.q.beta.Q2$s==s,]$c
        
        if(i!=state && j==end.state) {
          
          if((a==i && b==end.state)) {
            haz.ab <- 0
            for (u in 1:p) {
              haz.ab <- haz.ab + (data[mk,(u+2)]*beta[[u]][a, b])
            }
            haz.ab <- exp(haz.ab)
            
            d.q.beta.Q2[[c]][[s]][a,b] <- haz.ab * data[mk,(c+2)] * Q[state,i] 
          }
          
          if((a==state && b==i)) {
            d.q.beta.Q2[[c]][[s]][a,b] <- Q[i,end.state] * data[mk,(c+2)] * haz.ab 
          }
        }
        s <- s+1
      }
    }
    
    # 4
    d.beta.beta.Q2 <- list(list(list()))
    for (s in 1:p) {
      d.beta.beta.Q2[[s]] <- list(list())
    }
    for (s in 1:p) {
      for (r in 1:(nrow.Q*(nrow.Q-1))) {
        d.beta.beta.Q2[[s]][[r]] <- list()
      }
    }
    for (s in 1:p) {
      for (r in 1:(nrow.Q*(nrow.Q-1))) {
        for (u in 1:p) {
          d.beta.beta.Q2[[s]][[r]][[u]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
          diag(d.beta.beta.Q2[[s]][[r]][[u]]) <- NA
        }
      }
    }
    
    s <- 1
    ind.d.beta.beta.Q2 <- data.frame()
    for (c in 1:p) {
      for (r in 1:p) {
        for (a in 1:nrow.Q) {
          for (b in c(1:nrow.Q)[c(1:nrow.Q)!=a]) {
            
            ind.d.beta.beta.Q2 <- rbind(ind.d.beta.beta.Q2, c(s,a,b,c,r))
            s <- s+1
          }
        }
      }
    }
    colnames(ind.d.beta.beta.Q2) <- c('s','a','b','c','r')
    
    
    s <- 1
    for (i in 1:nrow.Q) {
      for (j in c(1:nrow.Q)[c(1:nrow.Q)!=i]) {
        a <- ind.d.beta.beta.Q2[ind.d.beta.beta.Q2$s==s,]$a
        b <- ind.d.beta.beta.Q2[ind.d.beta.beta.Q2$s==s,]$b
        c <- ind.d.beta.beta.Q2[ind.d.beta.beta.Q2$s==s,]$c
        r <- ind.d.beta.beta.Q2[ind.d.beta.beta.Q2$s==s,]$r
        
        if(i!=state && j==end.state) {
          
          if((a==i && b==end.state) ||  (a==state && b==i)) {
            d.beta.beta.Q2[[c]][[s]][[r]][a,b] <- Q[i,end.state] * data[mk,(r+2)] * data[mk,(c+2)] * Q[state,i]
          }
          
        }
        s <- s+1
      }
    }
    
    ####################################### 
    ####################################### Third situation
    # 1
    d.q.q.Q2 <- list()
    for (s in 1:(nrow.Q*(nrow.Q-1))) {
      d.q.q.Q2[[s]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
      diag(d.q.q.Q2[[s]]) <- NA
    }
    
    s <- 1
    ind.d.q.q.Q2 <- data.frame()
    for (a in 1:nrow.Q) {
      for (b in c(1:nrow.Q)[c(1:nrow.Q)!=a]) {
        ind.d.q.q.Q2 <- rbind(ind.d.q.q.Q2, c(s,a,b))
        s <- s+1
      }
    }
    colnames(ind.d.q.q.Q2) <- c('s','a','b')
    
    s <- 1
    for (i in 1:nrow.Q) {
      for (j in c(1:nrow.Q)[c(1:nrow.Q)!=i]) {
        a <- ind.d.q.q.Q2[ind.d.q.q.Q2$s==s,]$a
        b <- ind.d.q.q.Q2[ind.d.q.q.Q2$s==s,]$b
        
        if(i==state && j==end.state) {
          
          if(a==b && b==state){
            haz.state.state <- 0
            for (u in 1:p) {
              haz.state.state <- haz.state.state + (data[mk,(u+2)]*beta[[u]][state, state])
            }
            haz.state.state <- exp(haz.state.state)
            
            haz.state.end.state <- 0
            for (u in 1:p) {
              haz.state.end.state <- haz.state.end.state + (data[mk,(u+2)]*beta[[u]][state, end.state])
            }
            haz.state.end.state <- exp(haz.state.end.state)
            
            d.q.q.Q2[[s]][a,b] <- haz.state.state * haz.state.end.state
          }
          
          if(a==b && b==end.state){
            haz.end.state.end.state <- 0
            for (u in 1:p) {
              haz.end.state.end.state <- haz.end.state.end.state + (data[mk,(u+2)]*beta[[u]][end.state, end.state])
            }
            haz.end.state.end.state <- exp(haz.end.state.end.state)
            
            haz.state.end.state <- 0
            for (u in 1:p) {
              haz.state.end.state <- haz.state.end.state + (data[mk,(u+2)]*beta[[u]][state, end.state])
            }
            haz.state.end.state <- exp(haz.state.end.state)
            
            d.q.q.Q2[[s]][a,b] <- haz.end.state.end.state * haz.state.end.state
          }
        }
        s <- s+1
      }
    }
    
    # 2
    d.beta.q.Q2 <- list(list())
    for (s in 1:(nrow.Q*(nrow.Q-1))) {
      d.beta.q.Q2[[s]] <- list()
    }
    for (s in 1:(nrow.Q*(nrow.Q-1))) {
      for (r in 1:p) {
        d.beta.q.Q2[[s]][[r]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
        diag(d.beta.q.Q2[[s]][[r]]) <- NA
      }
    }
    
    s <- 1
    ind.d.beta.q.Q2 <- data.frame()
    for (c in 1:p) {
      for (a in 1:nrow.Q) {
        for (b in c(1:nrow.Q)[c(1:nrow.Q)!=a]) {
          ind.d.beta.q.Q2 <- rbind(ind.d.beta.q.Q2, c(s,a,b,c))
          s <- s+1
        }
      }
    }
    colnames(ind.d.beta.q.Q2) <- c('s','a','b','c')
    
    
    s <- 1
    for (i in 1:nrow.Q) {
      for (j in c(1:nrow.Q)[c(1:nrow.Q)!=i]) {
        a <- ind.d.beta.q.Q2[ind.d.beta.q.Q2$s==s,]$a
        b <- ind.d.beta.q.Q2[ind.d.beta.q.Q2$s==s,]$b
        c <- ind.d.beta.q.Q2[ind.d.beta.q.Q2$s==s,]$c
        
        if(i==state && j==end.state) {
          
          if(a==b && b==state) {
            haz.state.end.state <- 0
            for (u in 1:p) {
              haz.state.end.state <- haz.state.end.state + (data[mk,(u+2)]*beta[[u]][state, end.state])
            }
            haz.state.end.state <- exp(haz.state.end.state)
            
            d.beta.q.Q2[[s]][[c]][a,b] <- data[mk,(c+2)] * Q[state,state] * haz.state.end.state
          }
          
          if(a==b && b==end.state) {
            haz.state.end.state <- 0
            for (u in 1:p) {
              haz.state.end.state <- haz.state.end.state + (data[mk,(u+2)]*beta[[u]][state, end.state])
            }
            haz.state.end.state <- exp(haz.state.end.state)
            
            d.beta.q.Q2[[s]][[c]][a,b] <- data[mk,(c+2)] * Q[end.state,end.state] * haz.state.end.state
          }
          
          if(a==state && b==end.state) {
            haz.state.end.state <- 0
            for (u in 1:p) {
              haz.state.end.state <- haz.state.end.state + (data[mk,(u+2)]*beta[[u]][state, end.state])
            }
            haz.state.end.state <- exp(haz.state.end.state)
            
            d.beta.q.Q2[[s]][[c]][a,b] <- (Q[state,state] + Q[end.state,end.state]) * data[mk,(c+2)] * haz.state.end.state
          }
        }
        s <- s+1
      }
    }
    
    # 3
    d.q.beta.Q2 <- list(list())
    for (s in 1:p) {
      d.q.beta.Q2[[s]] <- list()
    }
    for (s in 1:p) {
      for (r in 1:(nrow.Q*(nrow.Q-1))) {
        d.q.beta.Q2[[s]][[r]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
        diag(d.q.beta.Q2[[s]][[r]]) <- NA
      }
    }
    
    s <- 1
    ind.d.q.beta.Q2 <- data.frame()
    for (a in 1:nrow.Q) {
      for (b in c(1:nrow.Q)[c(1:nrow.Q)!=a]) {
        for (c in 1:p) {
          ind.d.q.beta.Q2 <- rbind(ind.d.q.beta.Q2, c(s,a,b,c))
          s <- s+1
        }
      }
    }
    colnames(ind.d.q.beta.Q2) <- c('s','a','b','c')
    
    
    s <- 1
    for (i in 1:nrow.Q) {
      for (j in c(1:nrow.Q)[c(1:nrow.Q)!=i]) {
        a <- ind.d.q.beta.Q2[ind.d.q.beta.Q2$s==s,]$a
        b <- ind.d.q.beta.Q2[ind.d.q.beta.Q2$s==s,]$b
        c <- ind.d.q.beta.Q2[ind.d.q.beta.Q2$s==s,]$c
        
        if(i==state && j==end.state) {
          
          if((a==i && b==end.state)) {
            haz.ab <- 0
            for (u in 1:p) {
              haz.ab <- haz.ab + (data[mk,(u+2)]*beta[[u]][a, b])
            }
            haz.ab <- exp(haz.ab)
            
            d.q.beta.Q2[[c]][[s]][a,b] <- haz.ab * data[mk,(c+2)] * Q[state,i] 
          }
          
          if((a==state && b==i)) {
            d.q.beta.Q2[[c]][[s]][a,b] <- Q[i,end.state] * data[mk,(c+2)] * haz.ab 
          }
        }
        s <- s+1
      }
    }
    
    # 4
    d.beta.beta.Q2 <- list(list(list()))
    for (s in 1:p) {
      d.beta.beta.Q2[[s]] <- list(list())
    }
    for (s in 1:p) {
      for (r in 1:(nrow.Q*(nrow.Q-1))) {
        d.beta.beta.Q2[[s]][[r]] <- list()
      }
    }
    for (s in 1:p) {
      for (r in 1:(nrow.Q*(nrow.Q-1))) {
        for (u in 1:p) {
          d.beta.beta.Q2[[s]][[r]][[u]] <- matrix(0, nrow = nrow.Q, ncol = nrow.Q)
          diag(d.beta.beta.Q2[[s]][[r]][[u]]) <- NA
        }
      }
    }
    
    s <- 1
    ind.d.beta.beta.Q2 <- data.frame()
    for (c in 1:p) {
      for (r in 1:p) {
        for (a in 1:nrow.Q) {
          for (b in c(1:nrow.Q)[c(1:nrow.Q)!=a]) {
            
            ind.d.beta.beta.Q2 <- rbind(ind.d.beta.beta.Q2, c(s,a,b,c,r))
            s <- s+1
          }
        }
      }
    }
    colnames(ind.d.beta.beta.Q2) <- c('s','a','b','c','r')
    
    
    s <- 1
    for (i in 1:nrow.Q) {
      for (j in c(1:nrow.Q)[c(1:nrow.Q)!=i]) {
        a <- ind.d.beta.beta.Q2[ind.d.beta.beta.Q2$s==s,]$a
        b <- ind.d.beta.beta.Q2[ind.d.beta.beta.Q2$s==s,]$b
        c <- ind.d.beta.beta.Q2[ind.d.beta.beta.Q2$s==s,]$c
        r <- ind.d.beta.beta.Q2[ind.d.beta.beta.Q2$s==s,]$r
        
        if(i==state && j==end.state) {
          
          if(a==b && b==state) {
            haz.state.end.state <- 0
            for (u in 1:p) {
              haz.state.end.state <- haz.state.end.state + (data[mk,(u+2)]*beta[[u]][state, end.state])
            }
            haz.state.end.state <- exp(haz.state.end.state)
            
            d.beta.beta.Q2[[c]][[s]][[r]][a,b] <- Q[state,state] * Q[state,end.state] * data[mk,(r+2)] * data[mk,(c+2)] * haz.state.end.state
          }
          
          if(a==b && b==end.state) {
            haz.state.end.state <- 0
            for (u in 1:p) {
              haz.state.end.state <- haz.state.end.state + (data[mk,(u+2)]*beta[[u]][state, end.state])
            }
            haz.state.end.state <- exp(haz.state.end.state)
            
            d.beta.beta.Q2[[c]][[s]][[r]][a,b] <- Q[end.state,end.state] * Q[state,end.state] * data[mk,(r+2)] * data[mk,(c+2)] * haz.state.end.state
          }
          
          if(a==state && b==end.state) {
            haz.state.end.state <- 0
            for (u in 1:p) {
              haz.state.end.state <- haz.state.end.state + (data[mk,(u+2)]*beta[[u]][state, end.state])
            }
            haz.state.end.state <- exp(haz.state.end.state)
            
            d.beta.beta.Q2[[c]][[s]][[r]][a,b] <- 2*((Q[state,state] + Q[end.state,end.state]) * Q[state,end.state] * data[mk,(r+2)] * haz.state.end.state)
          }
          
        }
        s <- s+1
      }
    }
    
    
    vec.d.q.q.Q <- as.vector(t(do.call(rbind, d.q.q.Q)))
    vec.d.q.q.Q <- vec.d.q.q.Q[!is.na(vec.d.q.q.Q)]
    
    
    vec.d.beta.q.Q <- vector()
    for (s in 1:(nrow.Q*(nrow.Q-1))) {
      vec.d.beta.q.Q <- append(vec.d.beta.q.Q, t(do.call(rbind, d.beta.q.Q[[s]])))
    }
    vec.d.beta.q.Q <- vec.d.beta.q.Q[!is.na(vec.d.beta.q.Q)]
    
    
    vec.d.q.beta.Q <- vector()
    for (s in 1:p) {
      vec.d.q.beta.Q <- append(vec.d.q.beta.Q, t(do.call(rbind, d.q.beta.Q[[s]])))
    }
    vec.d.q.beta.Q <- vec.d.q.beta.Q[!is.na(vec.d.q.beta.Q)]
    
    
    vec.d.beta.beta.Q <- vector()
    for (s in 1:p) {
      for (r in 1:(nrow.Q*(nrow.Q-1))) {
        vec.d.beta.beta.Q <- append(vec.d.beta.beta.Q, t(do.call(rbind, d.beta.beta.Q[[s]][[r]])))
      }
    }
    vec.d.beta.beta.Q <- vec.d.beta.beta.Q[!is.na(vec.d.beta.beta.Q)]
    
    
    
    
    
    hessian.up <- cbind(matrix(vec.d.q.q.Q, nrow=nrow.Q*(nrow.Q-1), ncol=nrow.Q*(nrow.Q-1), byrow=TRUE), 
                        matrix(vec.d.beta.q.Q, nrow=nrow.Q*(nrow.Q-1), ncol=p*(nrow.Q*(nrow.Q-1)), byrow=TRUE))
    hessian.down <- cbind(matrix(vec.d.q.beta.Q, nrow=p*(nrow.Q*(nrow.Q-1)), ncol=nrow.Q*(nrow.Q-1), byrow=TRUE), 
                          matrix(vec.d.beta.beta.Q, nrow=p*(nrow.Q*(nrow.Q-1)), ncol=p*(nrow.Q*(nrow.Q-1)), byrow=TRUE))
    hessian.part.Q2 <- rbind(hessian.up, hessian.down)   
    
    hessian <- (hessian.part.Q*data[mk,2]) + (0.5*hessian.part.Q*(data[mk,2]^2))
    
    hes.list[[mk]] <- hessian
  }  
  sum.hessian <- Reduce('+', ll.foreach)
  stopCluster(cl)
  
  
  # fisher_info <- ginv(sum.hessian)
  # # if -2*log.lik
  # # fisher_info<-solve(0.5*res$hessian)
  # prop_sigma <- sqrt(diag(fisher_info))
  # prop_sigma <- diag(prop_sigma)
  # upper <- par+1.959964*prop_sigma
  # lower <- par-1.959964*prop_sigma
  # interval <- data.frame(value=par, upper=upper, lower=lower)
  
  return(sum.hessian)  
}

data <- readRDS('validate_data.RData')
par <- readRDS('opt_par.RData')
p <- ncol(data)-6
nrow.Q <- 8
hes <- hess.function(data=data, par=par )

saveRDS(hes, "hess.RData")
