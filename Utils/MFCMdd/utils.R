library(Rcpp)
library(Rcpp)
library(RcppArmadillo)




# unifInit
UInit <- function(data, n, k, m) {
  out <- matrix(0, n, k)
  H0_index <- sample.int(n=n, size = k, replace = FALSE)
  H <- data[H0_index,]
  H0 <- euclidean_distance_medoid(data, H, n, k, Square = FALSE) 
    
  for (i in 1:n) {
    d <- H0[i, ]
    if (min(d) == 0) {
      mm <- which.min(d)
      out[i, mm] <- 1
    } else {
      for (j in 1:k) {
        out[i, j] <- (1 / H0[i, j])^(1 / (m - 1)) / sum((1 / H0[i, ]))^(1 / (m - 1))
      }
    }
  }
  return(out)
}


# euclidean_distance_medoid
euclidean_distance_medoid <- function(data, H, n, k, Square = FALSE){
  n <- nrow(data)
  k <- nrow(H)
  out <- matrix(0, n, k)
  for (i in 1:n) {
    for (j in 1:k) {
      out[i, j] <- sum((data[i, ] - H[j, ])^2)
      if (Square) {
        out[i, j] <- sqrt(out[i, j])
      }
    }
  }
  
  return(out)
}



# memb_degree_medoid
memb_degree_medoid <- function(D, medoid, m, n, k) {
  n <- nrow(D)
  k <- ncol(D)
  out <- matrix(0, n, k)

  for (i in 1:n) {
    d <- D[i, ]
    index_temp_med <- i %in% medoid

    if (min(d) == 0) {
      mm <- which.min(d)
      out[i, mm] <- 1
    } else {
      for (j in 1:k) {
        out[i, j] <- (1 / D[i, j])^(1 / (m - 1)) / sum((1 / D[i, ]))^(1 / (m - 1))
      }
    }
  }

  return(out)
}




# Définition de la fonction mainFKM_med
mainFKM_med <- function(data, m, n, k, rs, conv, maxit, index, alpha, disp) {
  value <- numeric(rs)
  it <- numeric(rs)
  H <- matrix(0, nrow = k, ncol = n)
  D <- matrix(0, nrow = n, ncol = k)
  U <- matrix(0, nrow = n, ncol = k)
  func_opt <- 0
  
  H_opt <- matrix(0, nrow = k, ncol = n)
  U_opt <- matrix(0, nrow = n, ncol = k)
  U_old <- matrix(0, nrow = n, ncol = k)
  medoid_opt <- rep(0, k)
  min_med_const <- 10^5 * sum(data^2)
  min_med_old <- 0
  min_med <- 0
  
  medoid <- rep(0, k)
  index_temp_med <- TRUE
  value_temp_med <- TRUE
  
  ind <- 0
  ind_max <- 0
  nan_check <- 1
  
  for (r in 1:rs) {
    iter <- 0
    U <- UInit(data, n, k, m)
    U_old <- U
    prova <- TRUE
    
    while(prova && (iter < maxit))  {
      iter <- iter + 1
      U_old <- U
      
      dist_matrix <- as.matrix(data)
      for(c in 1:k) {
        min_med_i <- rep(Inf, n)
        
        for(i in 1:n) {
          if (all(is.finite(U[i,]))) {
            min_med_i[i] <- sum(U[,c]^m * dist_matrix[i,]^2)
          }
        }
        if (all(min_med_i == Inf)) {break}
        
        medoid[c] <- which.min(min_med_i)
        H[c,] <- data[medoid[c],]
      }
      
      D <- euclidean_distance_medoid(data, H, n, k)
      U <- memb_degree_medoid(D, medoid, m, n, k)

      if(disp) print(c(iter,sum(abs(U_old - U))))
      prova <- sum(abs(U_old - U)) > conv
    }
    
    if(nan_check == 0) {
      func <- NaN
    }else{
      func <- sum(U^m * D)
    }
    
    it[r] <- iter
    value[r] <- func
    if (is.finite(func)) {
      if (r == 1 || func < func_opt) {
        U_opt <- U
        H_opt <- H
        medoid_opt <- medoid
        func_opt <- func
      }
    }
  }
  
  # if (is.finite(func_opt)) {
  #   ind <- indices(index, data, U_opt, H_opt, m, n, k, p, exp(1), alpha)
  #   if(index == "PE" || index == "XB"){
  #     ind_max = -ind
  #   }else{
  #     ind_max = ind
  #   }
  # } else {
  #   ind <- NaN
  #   warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)")
  # }
  
  return(list(
    H = H_opt,
    U = U_opt,
    medoid = medoid_opt,
    iter = it,
    value = value,
    index = ind,
    index_max = ind_max,
    k = k
  ))
}



# Définition de la fonction mainFKM_med_U
mainFKM_med_U <- function(data, m, n, k, U, conv, maxit, index, alpha, disp) {
  iter <- 0
  H <- matrix(0, nrow = k, ncol = n)
  D <- matrix(0, nrow = n, ncol = k)
  U_old <- U
  prova <- TRUE
  
  value <- 0
  min_med_const <- 10^5 * sum(data^2)
  min_med_old <- 0
  min_med <- 0
  
  medoid <- rep(0, k)
  
  index_temp_med <- TRUE
  value_temp_med <- TRUE
  
  ind <- 0
  ind_max <- 0
  nan_check <- 1
  
  while(prova && (iter < maxit)) {

    iter <- iter + 1
    U_old <- U
    medoid <- rep(0, k)
    
    for (c in 1:k) {
      min_med_old <- min_med_const

      for (i in 1:n) {
        nan_check <- all(is.finite(U[i, ]))
        if (nan_check == 0) {
          break
        }
        min_med <- 0

        for (j in 1:n) {
          min_med <- min_med + sum(U[j, c]^m) * sum((data[i, ] - data[j, ])^2)
        }

        index_temp_med <- match(i, medoid, nomatch = 0) == 0
        value_temp_med <- min_med < min_med_old
        if (index_temp_med && value_temp_med) {
          min_med_old <- min_med
          medoid[c] <- i
        }
      }

      H[c, ] <- data[medoid[c], ]
    }
    
    D <- euclidean_distance_medoid(data, H, n, k)
    U <- memb_degree_medoid(D, medoid, m, n, k)
    
    if(disp) print(c(iter,U))
    prova <- sum(abs(U_old - U)) > conv
    
  }
  
  if (nan_check == 0) {
    ind <- NaN
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)")
  } else {
    value <- sum(U^m * D)
    ind <- indices(index, data, U, H, m, n, k, p, exp(1), alpha)
  }
  
  return(list(
    H = H,
    U = U,
    D = D,
    medoid = medoid + 1,
    iter = iter,
    value = value,
    index = ind,
    index_max = ind_max,
    k = k
  ))
}







## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



# Function index
indices <- function(type, X, U, H, m, n, k, p, b, alpha, distance = FALSE) {
  value <- 0
  
  if (type == "PC") {
    value <- partCoef(U, n)
  } else if (type == "PE") {
    value <- partEntropy(U, b, n)
  } else if (type == "MPC") {
    value <- partCoef_mod(U, n, k)
  } else if (type == "SIL") {
    value <- silhouette_internal(X, U, p, k, n, distance)
  } else if (type == "SIL.F") {
    value <- silhouetteFuzzy(X, U, alpha, p, k, n, distance)
  } else if (type == "XB") {
    value <- xie_beni(X, U, H, m, n, k)
  } else {
    stop("No match names.")
  }
  
  return(value)
}


# Fuzzy Silhouette Index :
silhouetteFuzzy <- function(X, U, alpha, p, k, n, distance = FALSE) {
  memb <- integer(n)
  count <- integer(k)
  D <- matrix(0, n, n)
  a <- numeric(n)
  b <- a
  sil_obj <- a
  B <- matrix(0, n, k)
  sil <- 0
  sil_F <- numeric(n)
  w <- numeric(n)
  ii <- integer(0)
  uu <- numeric(k - 1)
  u <- numeric(k)
  
  for (i in 1:n) {
    memb[i] <- which.max(U[i, ])
    count[memb[i]] <- count[memb[i]] + 1
  }
  
  if (!distance) {
    for (i in 1:(n - 1)) {
      for (i2 in (i + 1):n) {
        D[i, i2] <- sum((X[i, ] - X[i2, ])^2)
        D[i2, i] <- D[i, i2]
      }
    }
  } else {
    D <- X
  }
  
  for (i in 1:n) {
    for (j in 1:k) {
      for (i2 in 1:n) {
        if (memb[i2] == j) {
          B[i, j] <- B[i, j] + D[i, i2]
        }
      }
    }
  }
  
  for (i in 1:n) {
    for (j in 1:k) {
      if (memb[i] == j) {
        if (count[j] != 1) {
          B[i, j] <- B[i, j] / (count[j] - 1)
          a[i] <- B[i, j]
          B[i, j] <- max(B[i, ]) + 1
        }
      } else {
        B[i, j] <- B[i, j] / count[j]
      }
    }
    
    if (count[memb[i]] != 1) {
      b[i] <- min(B[i, ])
      sil_obj[i] <- (b[i] - a[i]) / max(a[i], b[i])
    }
  }
  
  for (i in 1:n) {
    u <- sort(U[i, ], decreasing = TRUE)
    uu <- u[2:k]
    w[i] <- (max(u) - max(uu))^alpha
  }
  
  sil <- sum(w * sil_obj) / sum(w)
  return(sil)
}



#Partition Coefficient :
partCoef <- function(U, n) {
  out <- sum(U^2) / n
  return(out)
}


#Partition Entropy Index :
partEntropy <- function(U, b, n) {
  out <- 0
  eps <- .Machine$double.eps
  logBase <- U
  
  U[U < eps] <- eps
  
  logBase <- log(U) / log(b)
  out <- -sum(U * logBase) / n
  
  return(out)
}

#Modified Partition Coefficient :
partCoef_mod <- function(U, n, k) {
  out <- 1 - k / (k - 1) * (1 - partCoef(U, n))
  return(out)
}

# euclidean_distance
euclidean_distance <- function(data, H, n, k, Square = FALSE) {
  out <- matrix(0, n, k)
  
  for (i in 1:n) {
    for (j in 1:k) {
      out[i, j] <- sum((data[i, ] - H[j, ])^2)
      if (Square) {
        out[i, j] <- sqrt(out[i, j])
      }
    }
  }
  return(out)
}


#Xie and Beni Index :
xie_beni <- function(X, U, H, m, n, k) {
  D <- euclidean_distance(X, H, n, k)
  distH <- 10^10 * sum(H^2)
  out <- 0
  
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      if (sum((H[i, ] - H[j, ])^2) < distH) {
        distH <- sum((H[i, ] - H[j, ])^2)
      }
    }
  }
  
  out <- sum(U^m * D) / (n * distH)
  return(out)
}


#Silhouette Index for Internal Use :
silhouette_internal <- function(X, U, p, k, n, distance = FALSE) {
  memb <- integer(n)
  count <- integer(k)
  
  D <- matrix(0, n, n)
  
  a <- numeric(n)
  b <- numeric(n)
  sil_obj <- numeric(n)
  B <- matrix(0, n, k)
  
  sil <- 0
  
  for (i in 1:n) {
    memb[i] <- which.max(U[i, ])
    count[memb[i]] <- count[memb[i]] + 1
  }
  
  if (!distance) {
    for (i in 1:(n-1)) {
      for (i2 in (i+1):n) {
        D[i, i2] <- sum((X[i, ] - X[i2, ])^2)
        D[i2, i] <- D[i, i2]
      }
    }
  } else {
    D <- X
  }
  
  for (i in 1:n) {
    for (j in 1:k) {
      for (i2 in 1:n) {
        if (memb[i2] == j) {
          B[i, j] <- B[i, j] + D[i, i2]
        }
      }
    }
  }
  
  for (i in 1:n) {
    for (j in 1:k) {
      if (memb[i] == j) {
        if (count[j] != 1) {
          B[i, j] <- B[i, j] / (count[j] - 1)
          a[i] <- B[i, j]
          B[i, j] <- max(B[i, ]) + 1
        }
      } else {
        B[i, j] <- B[i, j] / count[j]
      }
    }
    
    if (count[memb[i]] != 1) {
      b[i] <- min(B[i, ])
      sil_obj[i] <- (b[i] - a[i]) / max(a[i], b[i])
    }
  }
  
  sil <- mean(sil_obj)
  return(sil)
}



cl.memb <- function(U) {
  if (missing(U))
    stop("The membership degree matrix U must be given")
  
  if (is.null(U))
    stop("The membership degree matrix U is empty")
  
  U <- as.matrix(U)
  
  if (any(is.na(U)))
    stop("The membership degree matrix U must not contain NA values")
  
  if (!is.numeric(U))
    stop("The membership degree matrix U must be numeric")
  
  n <- nrow(U)
  info.U <- cbind(max.col(U), apply(U, 1, max))
  
  if (is.null(rownames(U)))
    rownames(info.U) <- paste("Obj", 1:n, sep=" ")
  else
    rownames(info.U) <- rownames(U)
  
  colnames(info.U) <- c("Cluster", "Membership degree")
  return(info.U)
}






## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

m_lambdaInit <- function(weight,p){
  lambda <- numeric(p)
  
  if(weight=='sum'){
    lambda <- rep((1/p),p)
  }else if(weight=='prod'){
    lambda <- rep(1,p)
  }
  return(lambda)
}


m_H_medoid <- function(data, medoid_index, p){
  out <- list()
  for(l in 1:p){
    out[[l]] <- t(data[[l]][medoid_index,])
  }
  return(out)
}


m_memb_degree_medoid <- function(data, n, k, p, m, weight, s, L, medoid_index){
  
  K <- k
  H <- mm_H_medoid(data, medoid_index, p)
  U <- matrix(0, n, K)
  
  for (i in 1:n) {
    for (k in 1:K) {
      numerator <- 0
      denominator <- 0
      
      for (h in 1:K) {
        sum_num <- 0
        sum_denom <- 0
        
        for (j in 1:p) {
          sum_num <- sum_num + (L[j]^s) * H[[j]][i, k]
          sum_denom <- sum_denom + (L[j]^s) * H[[j]][i, h]
        }
        numerator <- numerator + (sum_num / sum_denom)^(1 / (m - 1))
      }
      if(i == medoid_index[k]){
        U[i, k] <- 1
      }else{
        U[i, k] <- (1 / numerator)
      }
    }
  }
  return(U)
  return(out_sum)
}



m_memb_weight_lambda <-  function(data, U, m, n, k, p, weight, s, medoid_index){
  
  H <- m_H_medoid(data, medoid_index, p)
  lambda <- numeric(p)
  
  if(weight=="sum") {
    for(l in 1:p){
      HD_j <- H[[l]]
      
      v_ntor_dtor <- 0
      for(h in 1:p){
        HD_h <- H[[h]]
        
        dtor <- ntor <- 0
        for(j in 1:k){
          for(i in 1:n){
            ntor <- ntor + (U[i,j]^m) * HD_j[i,j]
            dtor <- dtor + (U[i,j]^m) * HD_h[i,j]
          }
        }
        v_ntor_dtor <- v_ntor_dtor + (ntor/dtor)^(1 / (s - 1))
      } 
      lambda[l] <- v_ntor_dtor^(-1)
    }
  }else if(weight=="prod"){
    for(l in 1:p){
      HD_j <- H[[l]]
      
      v_ntor <- 1
      for(h in 1:p){
        HD_h <- H[[h]]
        
        dtor <- ntor <- 0
        for(j in 1:k){
          for(i in 1:n){
            ntor <- ntor + (U[i,j]^m) * HD_h[i,j]
            dtor <- dtor + (U[i,j]^m) * HD_j[i,j]
          }
        }
        v_ntor <- v_ntor * ntor
      } 
      lambda[l] <- (v_ntor^(1/p)) / dtor
    }
  }
  return(lambda)
}


m_func <- function(data, U, m, n, k, p, s, L, medoid_index){
  
  H <- m_H_medoid(data, medoid_index, p)

  dtor <- 0
  for(h in 1:p){
    for(j in 1:k){
      for(i in 1:n){
        dtor <- dtor + (U[i,j]^m) * H[[h]][i,j] * (L[h]^s)
      }
    }
  } 
  return(dtor)
}


# m_medoid <- function(data, U, m, n, k, p, s, L){
#   medoid <-  rep(0, k)
#   
#   for(j in 1:k){
#     arg <- rep(0,n)
#     
#     for(h in 1:n){
#       arg_i <- 0
#       for(i in 1:n){
#         arg_l <- 0
#         for(l in 1:p){
#           data_p <- data[[l]]
#           arg_l <- arg_l + ((L[l]^s) * data_p[i,h])
#         }
#         arg_i <- arg_i + ((U[i,j]^m) * arg_l)
#       }
#       arg[h] <- arg_i
#     }
#     medoid[j] <- which.min(arg)  
#   }
#   return(medoid)
# }
m_medoid <- function(data, U, m, n, k, p, s, L){
  medoid <-  rep(0, k)
  
  for(j in 1:k){
    arg <- sapply(1:n, function(h) {
      arg_l <- sapply(1:p, function(l) {
        data_p <- data[[l]]
        (L[l]^s) * data_p[,h]
      })
      arg_i <- colSums((U[,j]^m) * arg_l)
      sum(arg_i)
    })
    medoid[j] <- which.min(arg)
  }
  return(medoid)
}



#  mainMFKM_med
mainMFKM_med <- function(data, m, n, p, k, conv, maxit, disp, weight, s) {

  iter <- min_med <- min_med_old <- 0
  nan_check <- 1
  
  H_opt <- list()
  U_opt <- U_old <- U <- matrix(0, nrow = n, ncol = k)
  min_med_const <- 10^5 * Reduce(`+`, lapply(data, function(df) sum(df^2)))
  
  medoid <- medoid_opt <- rep(0, k)
  index_temp_med <- value_temp_med <- TRUE
  
  # Initialization
  medoid_index <- sample.int(n=n, size = k, replace = FALSE)
  L <- m_lambdaInit(weight,p)
  U <- m_memb_degree_medoid(data, n, k, p, m, weight, s, L, medoid_index)
  Jold <- m_func(data, U, m, n, k, p, s, L, medoid_index)
  prova <- TRUE
    
  while(prova && (iter < maxit)){
    iter <- iter + 1
      
    medoid_index <- m_medoid(data, U, m, n, k, p, s, L)
    L <- m_memb_weight_lambda(data, U, m, n, k, p, weight, s, medoid_index)
    U <- m_memb_degree_medoid(data, n, k, p, m, weight, s, L, medoid_index)
    J <- m_func(data, U, m, n, k, p, s, L, medoid_index)
    
    if(disp) cat("iter: ",iter," J: ",J," epsi: ",abs(J-Jold),"\n")
    prova <- (abs(J-Jold)>conv)
    Jold <- J
  }
    
  func <- J
  H_opt <- m_H_medoid(data, medoid_index, p)
  U_opt <- U
  L_opt <- L
  
  
  return(list(
    H_opt = H_opt,
    U_opt = U_opt,
    medoid_opt = medoid_index,
    iter = iter,
    k = k,
    J = func,
    L_opt = L
  ))
}




# Définition de la fonction mainFKM_med_U
mainMFKM_med_U <- function(data, m, n, p, k, U, conv, maxit, index, alpha, disp) {
  iter <- 0
  H <- matrix(0, nrow = k, ncol = p)
  D <- matrix(0, nrow = n, ncol = k)
  U_old <- U
  prova <- TRUE
  
  value <- ind <- ind_max <- min_med_old <- min_med <- 0
  min_med_const <- 10^5 * sum(data^2)

  medoid <- rep(0, k)
  
  index_temp_med <- TRUE
  value_temp_med <- TRUE
  
  nan_check <- 1
  
  while (prova && (iter < maxit)) {
    iter <- iter + 1
    U_old <- U
    medoid <- rep(0, k)
    
    for (c in 1:k) {
      min_med_old <- min_med_const
      
      for (i in 1:n) {
        nan_check <- all(is.finite(U[i, ]))
        if (nan_check == 0) {
          break
        }
        min_med <- 0
        
        for (j in 1:n) {
          min_med <- min_med + sum(U[j, c]^m) * sum((data[i, ] - data[j, ])^2)
        }
        
        index_temp_med <- match(i, medoid, nomatch = 0) == 0
        value_temp_med <- min_med < min_med_old
        
        if (index_temp_med && value_temp_med) {
          min_med_old <- min_med
          medoid[c] <- i
        }
      }
      
      H[c, ] <- data[medoid[c], ]
    }
    
    D <- euclidean_distance_medoid(data, H, n, k)
    U <- memb_degree_medoid(D, medoid, m, n, k, p)
    
    if(disp) print(c(iter,U))
    prova <- sum(abs(U_old - U)) > conv
    
  }
  
  if (nan_check == 0) {
    ind <- NaN
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)")
  } else {
    value <- sum(U^m * D)
    ind <- indices(index, data, U, H, m, n, k, p, exp(1), alpha)
  }
  
  return(list(
    H = H,
    U = U,
    D = D,
    medoid = medoid + 1,
    iter = iter,
    value = value,
    index = ind,
    index_max = ind_max,
    k = k
  ))
}






## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

mm_lambdaInit <- function(weight, p, k){
  lambda <- matrix(0, nrow = k, ncol = p)
  for(j in 1:k){
    if(weight=='sum'){
      lambda[j,] <- rep((1/p),p)
    }else if(weight=='prod'){
      lambda[j,] <- rep(1,p)
    }
  }
  rownames(lambda) <- paste0("K_",1:k)
  colnames(lambda) <- paste0("P_",1:p)
  return(lambda)
}


mm_H_medoid <- function(data, medoid_index, p){
  out <- list()
  for(l in 1:p){
    out[[l]] <- t(data[[l]][medoid_index,])
  }
  return(out)
}


mm_memb_degree_medoid <- function(data, n, k, p, m, weight, s, L, medoid_index){
  
  K <- k
  H <- mm_H_medoid(data, medoid_index, p)
  U <- matrix(0, n, K)
  
  for (i in 1:n) {
    for (k in 1:K) {
      numerator <- 0
      denominator <- 0
      
      for (h in 1:K) {
        sum_num <- 0
        sum_denom <- 0
        
        for (j in 1:p) {
          sum_num <- sum_num + (L[k, j]^s) * H[[j]][i, k]
          sum_denom <- sum_denom + (L[h, j]^s) * H[[j]][i, h]
        }
        numerator <- numerator + (sum_num / sum_denom)^(1 / (m - 1))
      }
      if(i == medoid_index[k]){
        U[i, k] <- 1
      }else{
        U[i, k] <- (1 / numerator)
      }
    }
  }
  return(U)
}


mm_memb_weight_lambda <-  function(data, U, m, n, k, p, weight, s, medoid_index){

  H <- mm_H_medoid(data, medoid_index, p)
  lambda <- matrix(0, nrow = k, ncol = p)
  
  if(weight=="sum") {
    for(l in 1:p){
      HD_j <- H[[l]]
      for(j in 1:k){
        out_h <- 0
        for(h in 1:p){
          HD_h <- H[[h]]
          out_num <- out_deno <- 0
          for(i in 1:n){
            out_num <- out_num + (U[i,j]^m) * HD_j[i,j]
            out_deno <- out_deno + (U[i,j]^m) * HD_h[i,j]
          }
          out_h <-  out_h + (out_num/out_deno)^(1/(s-1))
        }
        lambda[j,l] <- out_h^(-1)
      }
    }
  }else if(weight=="prod"){
    for(l in 1:p){
      HD_j <- H[[l]]
      for(j in 1:k){
        
        out_h <- 1
        for(h in 1:p){
          HD_h <- H[[h]]
          out_num <- 0
          for(i in 1:n){
            out_num <-  out_num + U[i,j]^m * HD_h[i,j]
          }
          out_h <- out_h * out_num
        }
        out_deno <- 0
        for(i in 1:n){
          out_deno <-  out_deno + U[i,j]^m * HD_j[i,j]
        }
        lambda[j,l] <- (out_h^(1/p))/out_deno
      }
    }
  }
  return(lambda)
}


mm_func <- function(data, U, m, n, k, p, s, L, medoid_index){
  
  H <- mm_H_medoid(data, medoid_index, p)
  dtor <- 0
    for(j in 1:k){
      out_i <- 0
      for(i in 1:n){
        out_p <- 0
        for(l in 1:p){
          out_p <- out_p + ((L[j,l]^s) * H[[l]][i,j])
        }
        out_i <- out_i + ((U[i,j]^m) *  out_p)  
      }
      dtor <- dtor + out_i
  } 
  return(dtor)
}


mm_medoid <- function(data, U, m, n, k, p, s, L){
  medoid <-  rep(0, k)
  
  for(j in 1:k){
    arg <- sapply(1:n, function(h) {
      arg_l <- sapply(1:p, function(l) {
        data_p <- data[[l]]
        (L[j,l]^s) * data_p[,h]
      })
      arg_i <- colSums((U[,j]^m) * arg_l)
      sum(arg_i)
    })
    medoid[j] <- which.min(arg)
  }
  return(medoid)
}

unifInit <- function(n, d) {
  matrix(runif(n * d), nrow = n, ncol = d)
}




#  mainMFKM_med
mainMMFKM_med <- function(data, m, n, p, k, conv, maxit, disp, weight, s) {
  
  iter <- min_med <- min_med_old <- 0
  nan_check <- 1
  
  H_opt <- list()
  U_opt <- U_old <- U <- matrix(0, nrow = n, ncol = k)
  
  medoid <- medoid_opt <- rep(0, k)
  index_temp_med <- value_temp_med <- TRUE
  
  
  medoid_index <- sample.int(n=n, size = k, replace = FALSE)
  L <- mm_lambdaInit(weight,p,k)
  U <- mm_memb_degree_medoid(data, n, k, p, m, weight, s, L, medoid_index)
  Jold <- mm_func(data, U, m, n, k, p, s, L, medoid_index)
  prova <- TRUE
  
  while(prova && (iter < maxit)) {
    iter <- iter + 1
    
    medoid_index <- mm_medoid(data, U, m, n, k, p, s, L)
    L <- mm_memb_weight_lambda(data, U, m, n, k, p, weight, s, medoid_index)
    U <- mm_memb_degree_medoid(data, n, k, p, m, weight, s, L, medoid_index)
    J <- mm_func(data, U, m, n, k, p, s, L, medoid_index)

    if(disp) cat("iter: ",iter," J: ",J," epsi: ",abs(J-Jold),"\n")
    prova <- (abs(J-Jold)>conv)
    Jold <- J
  }
  
  func <- J
  H_opt <- mm_H_medoid(data, medoid_index, p)
  U_opt <- U
  L_opt <- L
  
  return(list(
    H_opt = H_opt,
    U_opt = U_opt,
    medoid_opt = medoid_index,
    iter = iter,
    k = k,
    J = func,
    L_opt = L
  ))
}




# Définition de la fonction mainFKM_med_U
mainMMFKM_med_U <- function(data, m, n, p, k, U, conv, maxit, index, alpha, disp) {
  iter <- 0
  H <- matrix(0, nrow = k, ncol = p)
  D <- matrix(0, nrow = n, ncol = k)
  U_old <- U
  prova <- TRUE
  
  value <- ind <- ind_max <- min_med_old <- min_med <- 0
  min_med_const <- 10^5 * sum(data^2)
  
  medoid <- rep(0, k)
  
  index_temp_med <- TRUE
  value_temp_med <- TRUE
  
  nan_check <- 1
  
  while (prova && (iter < maxit)) {
    iter <- iter + 1
    U_old <- U
    medoid <- rep(0, k)
    
    for (c in 1:k) {
      min_med_old <- min_med_const
      
      for (i in 1:n) {
        nan_check <- all(is.finite(U[i, ]))
        if (nan_check == 0) {
          break
        }
        min_med <- 0
        
        for (j in 1:n) {
          min_med <- min_med + sum(U[j, c]^m) * sum((data[i, ] - data[j, ])^2)
        }
        
        index_temp_med <- match(i, medoid, nomatch = 0) == 0
        value_temp_med <- min_med < min_med_old
        
        if (index_temp_med && value_temp_med) {
          min_med_old <- min_med
          medoid[c] <- i
        }
      }
      
      H[c, ] <- data[medoid[c], ]
    }
    
    D <- euclidean_distance_medoid(data, H, n, k)
    U <- memb_degree_medoid(D, medoid, m, n, k, p)
    
    if(disp) print(c(iter,U))
    prova <- sum(abs(U_old - U)) > conv
    
  }
  
  if (nan_check == 0) {
    ind <- NaN
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)")
  } else {
    value <- sum(U^m * D)
    ind <- indices(index, data, U, H, m, n, k, p, exp(1), alpha)
  }
  
  return(list(
    H = H,
    U = U,
    D = D,
    medoid = medoid + 1,
    iter = iter,
    value = value,
    index = ind,
    index_max = ind_max,
    k = k
  ))
}





## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Calcul du Coefficient de Partition (PC)
PC2 <- function(U, m) {
  N <- length(U)
  return((1/N) * sum((1/m) * sum(U^m, na.rm = TRUE), na.rm = TRUE))
}

# Calcul de l'Entropie de Partition (PE)
PE2 <- function(U) {
  N <- length(U)
  return(-1/N * sum(U * log(U), na.rm = TRUE))
}

# Calcul du Coefficient de Partition Modifié (MPC)
MPC2 <- function(U, m) {
  c <- length(U)
  return(PC(U, m) * c)
}




