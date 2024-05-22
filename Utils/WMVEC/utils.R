Centroids_Initialization <- function(X, K) {
  centroids <- matrix(NA, nrow = 0, ncol = ncol(X)) # Initialize centroids
  
  for (i in 1:K) {
    if (i == 1) {
      d <- numeric() # Initialize distances vector
      
      for (j in 1:nrow(X)) {
        center <- X[j,]
        dis <- 0
        
        for (z in 1:nrow(X)) {
          # Calculate the distance except for the centroid
          if (z != j) {
            dis <- dis + norm(center - X[z,], type = "2")
          }
        }
        d <- c(d, dis)
      }
      
      idx <- which.min(d)
      centroids <- rbind(centroids, X[idx,])
      X <- X[-idx,]
    } else {
      if (nrow(centroids) == 1) {
        for (j in 1:nrow(centroids)) {
          center <- centroids[j,]
          dis <- rowSums((X - center)^2)
          idx <- which.max(dis)
          centroids <- rbind(centroids, X[idx,])
          X <- X[-idx,]
        }
      } else {
        Dis <- numeric()
        Idx <- numeric()
        
        for (j in 1:nrow(centroids)) {
          center <- centroids[j,]
          dis <- rowSums((X - center)^2)
          Dis <- c(Dis, sum(dis))
          idx <- which.max(dis)
          Idx <- c(Idx, idx)
        }
        
        idx <- which.min(Dis)
        centroids <- rbind(centroids, X[Idx[idx],])
        X <- X[-Idx[idx],]
      }
    }
  }
  
  return(centroids)
}





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
get_distance <- function(mode, view, data, nbFoc, Aj, F_update, alpha, 
                         beta, delta, features, R_or_M){
  
  Dis <- vector("list", view)
  if(mode == 1) {
    R <- R_or_M
  }
  
  if(mode == 2) {
    M <- R_or_M
  }
  
  for (i in 1:view) {
    Row <- nrow(data[[i]])
    ROW <- nrow(F_update[[i]])
    dis_temp <- matrix(0, nrow = Row, ncol = ROW)
    #dis_cell <- vector("list", features[i])
    dis_cell <- list()
    
    
    for (p in 1:features[i]) {
      dis_cell[[p]] <- matrix(0, nrow = Row, ncol = nbFoc)
    }
    row <- nrow(Aj[[i]])
    
    for(j in 2:row) {
      temp <- (data[[i]] - Aj[[i]][j,])^2

      if(sum(F_update[[i]][j,]) > 0) {
        card <- sum(F_update[[i]][j,])
      }else {
        card <- 0
      }
      
      if(mode == 1) {
        temp <- rowSums(temp)
        temp[temp == 0] <- 1e-10
        dis_temp[,j] <- (temp * R[i] * card^alpha)
      }
      
      if(mode == 2) {
        temp <- rowSums(temp)
        temp[temp == 0] <- 1e-10
        dis_temp[,j] <- (temp * (M[,j]^beta) * card^alpha)
      }
    }
    
    if(mode == 1) {
      dis_temp[,1] <- delta[i]^2 * R[i]
      Dis[[i]] <- dis_temp
    }
    
    if(mode == 2) {
      dis_temp[,1] <- delta[i]^2 * (M[,1]^beta)
      Dis[[i]] <- dis_temp
    }
  }
  
  return(Dis)
}




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

update_Aj <- function(view, cluster, features, Aj, F, center, nbFoc) {
  for (i in 1:view) {
    new_center <- matrix(0, nrow = nbFoc, ncol = features[i])
    
    for(j in 1:nrow(F)) {
      
      if(sum(F[j,]) != 0) {
        temp1 <- 0
        
        for(k in 1:cluster) {
          temp1 <- temp1 + F[j,k] * center[[i]][k,]
        }
        
        new_center[j,] <- temp1 / sum(F[j,])
      }
    }
    Aj[[i]] <- new_center
  }
  return(list(Aj, new_center))
}





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

update_M <- function(view, data, dis, beta, nbFoc) {
  
  # Update the belief matrix M
  row <- nrow(data[[1]])
  M <- matrix(0, nrow = row, ncol = nbFoc)
  index <- 1 / (beta - 1)
  D <- matrix(0, nrow = row, ncol = nbFoc)
  
  for (i in 1:view) {
    D <- D + dis[[i]]
  }
  
  # Update M
  m <- matrix(0, nrow = row, ncol = nbFoc - 1)
  
  for(i in 1:row) {
    vect0 <- D[i, 2:ncol(D)]
    
    for(j in 2:nbFoc) {
      if(sum(D[, j]) != 0){
        vect1 <- ((D[i, j] * rep(1, nbFoc - 1)) / vect0)^(1 / (beta - 1))
        
        vect1[is.infinite(vect1)] <- 0
        vect3 <- vect1
        
        m[i, j-1] <- 1 / (sum(vect3) + (D[i, j] / D[i, 1])^(1 / (beta - 1))) ####
      }
    }
  }
  m <- cbind(1 - rowSums(m), m)
  M <- m
  M[M[,1] < 0] <- 0
  return(M)
}




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

update_R <- function(view, dis, lambda) {
  R <- rep(1/view, view)
  temp <- c()
  
  for (i in 1:view) {
    F_i <- sum(dis[[i]])
    temp <- c(temp, exp(-F_i/lambda))
  }
  
  
  for (i in 1:view) {
    R[i] <- temp[i] / sum(temp)
  }
  
  return(R)
}




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
update_V <- function(view, cluster, alpha, beta, data, M, R, 
                     F_update, features) {
  
  V <- list()
  
  for(i in 1:view) {
    v_i <- matrix(0, nrow = cluster, ncol = features[i])
    
    for(p in 1:features[i]) {
      B <- matrix(0, nrow = cluster, ncol = 1)
      
      for(j in 1:cluster) {
        pos <- c()
        
        for(k in 1:nrow(F_update[[i]])) {
          if(F_update[[i]][k,j] == 1) {
            pos <- c(pos, k)
          }
        }
        
        B_X <- 0
        
        for (n in 1:length(pos)) {
          card <- sum(F_update[[i]][pos[n],])
          r <- R[i]
          aj_dis <- card^(alpha-1) * r * data[[i]][,p] * (M[,pos[n]]^beta)
          B_X <- B_X + sum(aj_dis)
        }
        B[j] <- B_X
      }
      
      H <- matrix(0, nrow = cluster, ncol = cluster)
      
      for (c in 1:cluster) {
        for (k in 1:cluster) {
          loc <- c()
          for (n in 1:nrow(F_update[[i]])) {
            if (F_update[[i]][n,c] == 1 && F_update[[i]][n,k] == 1) {
              loc <- c(loc, n)
            }
          }
          
          H_ck <- 0
          
          for (n in 1:length(loc)) {
            card <- sum(F_update[[i]][loc[n],])
            r <- R[i]
            aj_all_dis <- card^(alpha-2) * r * M[,loc[n]]^beta
            H_ck <- H_ck + sum(aj_all_dis)
          }
          H[c,k] <- H_ck
        }
      }
      v <- solve(H) %*% B
      v_i[,p] <- v
    }
    V[[i]] <- v_i
  }
  return(V)
}




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
update_jaccard <- function(view, lambda, R, dis) {
  part_1 <- 0
  
  for (i in 1:view) {
    part_1 <- part_1 + R[i] * sum(dis[[i]])
  }
  
  part_2 <- 0
  
  for (i in 1:length(R)) {
    part_2 <- part_2 + lambda * sum((R[i] + 1e-4) * log(R[i] + 1e-3))
  }
  
  jaccard <- part_1 + part_2
  return(jaccard)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
extractMass <- function(mass,F,g=NULL,S=NULL,method,crit=NULL,Kmat=NULL,trace=NULL,
                        D=NULL,W=NULL,J=NULL,param=NULL){
  
  n<-nrow(mass)
  c<-ncol(F)
  if(any(F[1,]==1)){
    F<-rbind(rep(0,c),F) 
    mass<-cbind(rep(0,n),mass)
  }
  f<-nrow(F)
  card<-rowSums(F)
  
  conf<-mass[,1]            
  C<- 1/(1-conf)
  mass.n<- C*mass[,2:f]   
  pl<- mass%*% F          
  pl.n<- C*pl             
  p <-pl/rowSums(pl)      
  bel<- mass[,card==1]   
  bel.n<-C*bel            
  y.pl<-max.col(pl)       
  y.bel<-max.col(bel)     
  Y<-F[max.col(mass),]    
  # non dominated elements
  Ynd<-matrix(0,n,c)
  for (i in 1:n){
    ii<-which(pl[i,]>= bel[i,y.bel[i]])
    Ynd[i,ii]<-1
  }
  
  P<-F/card
  P[1,]<- rep(0,c)
  betp<- mass %*% P   
  betp.n<- C* betp    
  
  lower.approx<-vector(mode='list',length=c)
  upper.approx<-vector(mode='list',length=c)
  lower.approx.nd<-vector(mode='list',length=c)
  upper.approx.nd<-vector(mode='list',length=c)
  nclus<-rowSums(Y)
  outlier<-which(nclus==0) # outliers
  nclus.nd<-rowSums(Ynd)
  for(i in 1:c){
    upper.approx[[i]]<- which(Y[,i]==1)                 
    lower.approx[[i]]<- which((Y[,i]==1) & (nclus==1))  
    upper.approx.nd[[i]]<- which(Ynd[,i]==1)                 
    lower.approx.nd[[i]]<- which((Ynd[,i]==1) & (nclus.nd==1)) 
  }
  # Nonspecificity
  card<-c(c,card[2:f])
  Card<-matrix(card,n,f,byrow=TRUE)
  N <- sum(log(Card)*mass)/log(c)/n
  clus<-list(conf=conf,F=F,mass=mass,mass.n=mass.n,pl=pl,pl.n=pl.n,bel=bel,bel.n=bel.n,y.pl=y.pl,
             y.bel=y.bel,Y=Y,betp=betp,betp.n=betp.n,p=p,upper.approx=upper.approx,
             lower.approx=lower.approx,Ynd=Ynd,upper.approx.nd=upper.approx.nd,
             lower.approx.nd=lower.approx.nd,N=N,outlier=outlier,g=g,S=S,
             crit=crit,Kmat=Kmat,trace=trace,D=D,method=method,W=W,J=J,param=param)
  class(clus)<-"credpart"
  return(clus)
}


















###############################################################################
###############################################################################
###############################################################################

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
valid_external <- function(index1, c2) {
  
  N <- ncol(index1)
  Q <- sum(as.integer(c2) - c2) != 0
  Q <- Q || sum(sum(as.integer(index1) - index1)) != 0
  
  if (Q) {
    return(NULL)
  }
  
  Outs <- c()
  
  for (i in 1:N) {
    c1 <- index1[, i]
    C <- Contingency(c1, c2)
    n <- length(c1)
    nis <- sum(rowSums(C)^2)
    njs <- sum(colSums(C)^2)
    ns <- n * (n - 1) / 2
    sumC <- sum(C^2)
    sumij <- nis + njs
    R <- ns + sumC - sumij * 0.5
    nc <- (n * (n^2 + 1) - (n + 1) * nis - (n + 1) * njs + 2 * (nis * njs) / n) / (2 * (n - 1))
    
    if (ns == nc) {
      AR <- 0
    } else {
      AR <- (R - nc) / (ns - nc) # adjusted Rand - Hubert & Arabie 1985
    }
    
    Rand <- R / ns
    Jac <- (sumC - n) / (sumij - sumC - n) # Jaccard - Jain and Dubes 1988
    ni <- rowSums(C)
    ni <- ni * (ni - 1) / 2
    nis <- sum(ni)
    nj <- colSums(C)
    nj <- nj * (nj - 1) / 2
    njs <- sum(nj)
    FM <- 0.5 * (sumC - n) / sqrt(nis * njs) # FM - Fowlkes and Mallows 1983
    Outs <- c(Outs, c(Rand, AR, Jac, FM))
  }
  return(Outs)
}

Contingency <- function(Mem1, Mem2) {
  Cont <- matrix(0, nrow = max(Mem1), ncol = max(Mem2))
  for (i in 1:length(Mem1)) {
    Cont[Mem1[i], Mem2[i]] <- Cont[Mem1[i], Mem2[i]] + 1
  }
  return(Cont)
}




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
hmreduce <- function(A, CH, RH, LC, LR, SLC, SLR) {
  n <- nrow(A)
  
  coveredRows <- LR == 0
  coveredCols <- LC != 0
  
  r <- which(!coveredRows)
  c <- which(!coveredCols)
  m <- min(A[r, c])
  A[r, c] <- A[r, c] - m
  
  for (j in c) {
    for (i in SLC) {
      if (A[i, j] == 0) {
        if (RH[i] == 0) {
          RH[i] <- RH[n + 1]
          RH[n + 1] <- i
          CH[i] <- j
        }
        row <- A[i, ]
        colsInList <- -row[row < 0]
        if (length(colsInList) == 0) {
          l <- n + 1
        } else {
          l <- colsInList[row[colsInList] == 0]
        }
        A[i, l] <- -j
      }
    }
  }
  r <- which(coveredRows)
  c <- which(coveredCols)
  
  i_j <- which(A[r, c] <= 0, arr.ind = TRUE)
  i <- r[i_j[, 1]]
  j <- c[i_j[, 2]]
  
  for (k in 1:length(i)) {
    lj <- which(A[i[k], ] == -j[k])
    A[i[k], lj] <- A[i[k], j[k]]
    A[i[k], j[k]] <- 0
  }
  A[r, c] <- A[r, c] + m
  
  return(list(A = A, CH = CH, RH = RH))
}


hmflip <- function(A, C, LC, LR, U, l, r) {
  n <- nrow(A)
  while (TRUE) {
    C[l] <- r
    m <- which(A[r, ] == -l)
    
    A[r, m] <- A[r, l]
    
    A[r, l] <- 0
    if (LR[r] < 0) {
      U[n + 1] <- U[r]
      U[r] <- 0
      return(list(A = A, C = C, U = U))
    } else {
      l <- LR[r]
      A[r, l] <- A[r, n + 1]
      A[r, n + 1] <- -l
      r <- LC[l]
    }
  }
}


hminiass <- function(A) {
  n <- nrow(A)
  C <- rep(0, n)
  U <- rep(0, n + 1)
  LZ <- rep(0, n)
  NZ <- rep(0, n)
  
  for (i in 1:n) {
    lj <- n + 1
    j <- -A[i, lj]
    while (C[j] != 0) {
      lj <- j
      j <- -A[i, lj]
      if (j == 0) {
        break
      }
    }
    
    if (j != 0) {
      C[j] <- i
      A[i, lj] <- A[i, j]
      NZ[i] <- -A[i, j]
      LZ[i] <- lj
      A[i, j] <- 0
    } else {
      lj <- n + 1
      j <- -A[i, lj]
      while (j != 0) {
        r <- C[j]
        lm <- LZ[r]
        m <- NZ[r]
        while (m != 0) {
          if (C[m] == 0) {
            break
          }
          lm <- m
          m <- -A[r, lm]
        }
        
        if (m == 0) {
          lj <- j
          j <- -A[i, lj]
        } else {
          A[r, lm] <- -j
          A[r, j] <- A[r, m]
          NZ[r] <- -A[r, m]
          LZ[r] <- j
          A[r, m] <- 0
          C[m] <- r
          A[i, lj] <- A[i, j]
          NZ[i] <- -A[i, j]
          LZ[i] <- lj
          A[i, j] <- 0
          C[j] <- i
          
          break
        }
      }
    }
  }
  r <- rep(0, n)
  rows <- C[C != 0]
  r[rows] <- rows
  empty <- which(r == 0)
  
  U <- rep(0, n + 1)
  U[c(n + 1, empty)] <- c(empty, 0)
  
  return(list(A = A, C = C, U = U))
}


hminired <- function(A) {
  m <- nrow(A)
  n <- ncol(A)
  
  colMin <- apply(A, 2, min)
  A <- sweep(A, 2, colMin)
  
  rowMin <- apply(A, 1, min)
  A <- sweep(A, 1, rowMin, "-")
  
  i_j <- which(A == 0, arr.ind = TRUE)
  i <- i_j[, 1]
  j <- i_j[, 2]
  
  A <- cbind(A, rep(0, m))
  for (k in 1:n) {
    cols <- j[k == i]
    A[k, c(n + 1, cols)] <- c(-cols, 0)
  }
  
  return(A)
}


hungarian <- function(A) {
  m <- nrow(A)
  n <- ncol(A)
  
  if (m != n) {
    stop("HUNGARIAN: Cost matrix must be square!")
  }
  orig <- A
  A <- hminired(A)
  res <- hminiass(A)
  A <- res$A
  C <- res$C
  U <- res$U
  while (U[n + 1]) {
    LR <- rep(0, n)
    LC <- rep(0, n)
    CH <- rep(0, n)
    RH <- c(rep(0, n), -1)
    SLC <- c()
    r <- U[n + 1]
    LR[r] <- -1
    SLR <- r
    while (TRUE) {
      if (A[r, n + 1] != 0) {
        l <- -A[r, n + 1]
        if (A[r, l] != 0 & RH[r] == 0) {
          RH[r] <- RH[n + 1]
          RH[n + 1] <- r
          CH[r] <- -A[r, l]
        }
      } else {
        if (RH[n + 1] <= 0) {
          res <- hmreduce(A, CH, RH, LC, LR, SLC, SLR)
          A <- res$A
          CH <- res$CH
          RH <- res$RH
        }
        
        r <- RH[n + 1]
        l <- CH[r]
        CH[r] <- -A[r, l]
        if (A[r, l] == 0) {
          RH[n + 1] <- RH[r]
          RH[r] <- 0
        }
      }
      
      while (LC[l] != 0) {
        if (RH[r] == 0) {
          if (RH[n + 1] <= 0) {
            res <- hmreduce(A, CH, RH, LC, LR, SLC, SLR)
            A <- res$A
            CH <- res$CH
            RH <- res$RH
          }
          r <- RH[n + 1]
        }
        l <- CH[r]
        CH[r] <- -A[r, l]
        if(A[r, l] == 0) {
          RH[n + 1] <- RH[r]
          RH[r] <- 0
        }
      }
      if (C[l] == 0) {
        res <- hmflip(A, C, LC, LR, U, l, r)
        A <- res$A
        C <- res$C
        U <- res$U
        break
      } else {
        LC[l] <- r
        SLC <- c(SLC, l)
        r <- C[l]
        LR[r] <- l
        SLR <- c(SLR, r)
      }
    }
  }
  
  T <- sum(orig[as.logical(sparseMatrix(C, 1:ncol(orig), 
                                        x = rep(1, ncol(orig))))])
  
  return(list(C = C, T = T))
}


bestMap <- function(L1, L2) {
  L1 <- as.vector(L1)
  L2 <- as.vector(L2)
  if (length(L1) != length(L2)) {
    stop("size(L1) must size(L2)")
  }
  L1 <- L1 - min(L1) + 1 
  L2 <- L2 - min(L2) + 1
  
  nClass <- max(c(max(L1), max(L2)))
  G <- matrix(0, nrow = nClass, ncol = nClass)
  for (i in 1:nClass) {
    for (j in 1:nClass) {
      G[i, j] <- length(which(L1 == i & L2 == j))
    }
  }
  
  list_C_T <- hungarian(-G)
  c <- list_C_T[[1]]
  newL2 <- rep(0, nClass)
  for (i in 1:nClass) {
    newL2[L2 == i] <- c[i]
  }
  
  return(list(newL2 = newL2, c = c))
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
CalcMeasures <- function(Y, predY) {
  # result = [ACC,error_cnt]
  if (ncol(Y) != 1) {
    Y <- t(Y)
  }
  
  if (ncol(predY) != 1) {
    predY <- t(predY)
  }
  
  # bestMap
  predY <- bestMap(Y, predY)
  if (length(Y) != length(predY)) {
    predY <- t(predY)
  }
  
  error_cnt <- sum(Y != predY)
  AC <- length(which(Y == predY)) / length(Y)
  
  result <- c(AC, error_cnt)
  
  return(result)
}