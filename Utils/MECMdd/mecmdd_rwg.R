mecmdd.rwg<-function(Xlist, c, type='full', alpha=1, beta=1.5, delta=9, epsi=1e-4, 
                    disp=TRUE, gamma = 0.5, eta = 1, weight='sum', s=NULL){
  
  
  
  if(missing(Xlist) || is.null(Xlist) || !is.list(Xlist) || length(Xlist) <= 1){
    stop("The data set (Xlist) must be given, must not be empty, must be a list, and must contain more than one matrix")
  }
  
  tailles <- sapply(Xlist, dim)
  if(!all(tailles == tailles[,1])) {
    stop("Toutes les matrices dans Xlist n'ont pas la même taille")
  }
  
  if(weight=='prod'){s <- 1
  }else if(weight=='sum'){
    if(is.null(s)){
      s <- 2
    }
  }
  
  pmatrix <- length(Xlist)
  if(is.null(delta)){
    delta2 <- lapply(Xlist, function(D) quantile(D[upper.tri(D) | lower.tri(D)], 0.95))
    delta2 <- sum(unlist(delta2))
  }else if(length(delta) != pmatrix){
    delta2<-rep(delta^2,pmatrix)
    delta2 <- sum(delta2)
  }else{
    delta2<- delta^2
    delta2 <- sum(delta2)
  }
  
  F<-makeF(c=c,type=type)
  nF <- nrow(F)
  card<- rowSums(F[2:nF,])
  
  
  n <- nrow(Xlist[[1]])
  medoids <- sample.int(n, c, replace = FALSE)
  prototype<- numeric(nF-1)
  lambda <- lambdaInit_global(weight, pmatrix)
  
  pasfini<-TRUE
  Jold <- Inf
  iter<-0
  
  while(pasfini){
    iter<-iter+1
    
    # Imprecis
    for(j in 1:(nF-1)) {
      fj <- F[j+1, ]
      medoidsj <- medoids[which(fj != 0)]
      l_i <- numeric(n)
      
      if(sum(fj) == 1){
        prototype[j] <- medoidsj
      }else{
        for(i in 1:n){
          var_ij <-  (1/card[j]) * sum(sapply(seq_along(medoidsj), function(k)
            (sum(sapply(seq_along(Xlist), function(l) Xlist[[l]][i,medoidsj[k]])) -
               ((1/card[j]) * sum(sapply(seq_along(medoidsj), function(kk)
                 sum(sapply(seq_along(Xlist), function(l) Xlist[[l]][i,medoidsj[kk]]))))))^2))
          
          l_i[i] <- var_ij + eta * ((1/card[j]) * sum(sapply(seq_along(medoidsj), function(kk)
            sum(sapply(seq_along(Xlist), function(l) Xlist[[l]][i,medoidsj[kk]])))))  
        }
        prototype[j] <- which.min(l_i)
      }
    }
    
    # masse
    m <- matrix(0,n,nF-1)
    for(j in 1:(nF-1)){
      fj <- F[j+1, ]
      medoidsj <- medoids[which(fj != 0)]
      
      for(i in 1:n){
        num <- (card[j]^alpha * sum(sapply(seq_along(Xlist), 
                                           function(l) Xlist[[l]][i,prototype[j]] * (lambda[l])^s)))^(-1/(beta-1))
        den <- ((sum(sapply(seq(nF-1), function(k) 
          (card[k]^alpha * sum(sapply(seq_along(Xlist), 
            function(l) Xlist[[l]][i,prototype[k]] * (lambda[l])^s)))^(-1/(beta-1))  ))) + delta2^(-1/(beta-1)))
        m[i,j] <- num/den
        if (is.nan(m[i,j]) || is.infinite(m[i,j])) m[i,j] <- 1
      }
    }
    
    
    
    # weight
    if(weight=='sum'){
      for(l in 1:pmatrix){
        num <- sum(sapply(seq(n), function(i)
          sum(sapply(seq(nF-1), function(j)
            s * (card[j]^alpha) * m[i,j]^beta *  Xlist[[l]][i,prototype[j]]))))^(-1/(s-1))
        
        den <- sum(sapply(seq_along(Xlist), function(h)
          sum(sapply(seq(n), function(i)
            sum(sapply(seq(nF-1), function(j)
              s * (card[j]^alpha) * m[i,j]^beta *  Xlist[[h]][i,prototype[j]]))))^(-1/(s-1))))
        
        lambda[l] <- num/den
      }
    }
    
    if(weight=='prod'){
      for(l in 1:pmatrix){
        num <- prod(sapply(seq_along(Xlist), function(h)
          sum(sapply(seq(n), function(i)
            sum(sapply(seq(nF-1), function(j)
              (card[j]^alpha) * m[i,j]^beta *  Xlist[[h]][i,prototype[j]]))))^(1/pmatrix))
          ,na.rm = TRUE)
        
        den <- sum(sapply(seq(n), function(i)
          sum(sapply(seq(nF-1), function(j)
            (card[j]^alpha) * m[i,j]^beta *  Xlist[[l]][i,prototype[j]]))))
        
        lambda[l] <- num/den
      }
    }

    

    
    # medoids
    V <- matrix(0,n,nF-1)
    for(k in 1:(nF-1)) {
      fk <- F[k+1, ]
      if(sum(fk) == 1) {
        mk_beta <- m[,k]^beta
        for(i in 1:n){
          V[i,k] <- sum(rowSums(sapply(seq_along(Xlist), 
                                       function(l) Xlist[[l]][i,] * (lambda[l])^s))*mk_beta)
        }
      }
    }
    V_filtered <- V[, colSums(V) > 0]
    medoids <- sapply(1:ncol(V_filtered), function(i) which.min(V_filtered[,i]))
    
    
    mvide <- 1-rowSums(m)
    J <- sum(sapply(seq(n), function(i) 
      sum(sapply(seq(nF-1), function(j) card[j]^alpha * m[i,j]^beta * 
            sum(sapply(seq_along(Xlist), 
                function(l) Xlist[[l]][i,prototype[j]] * (lambda[l])^s)))))) +
      delta2*sum(mvide^beta, na.rm = TRUE)
    
    
    if(disp) print(c(iter,J))
    pasfini <- (abs(J-Jold)>epsi)
    Jold <- J
  }
  
  g <- g <- medoids
  m <- cbind(1-rowSums(m),m)
  clus<-extractMass(m,F,g=g,method="mecmdd",crit=j,
                    param=list(alpha=alpha,beta=beta,delta=delta,lambda=lambda))
  return(clus)
}
