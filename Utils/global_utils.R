weighted_silhouette <- function(labels, distance_matrices, weight_matrix, local) {
  n <- length(labels)
  num_matrices <- length(distance_matrices)
  
  if(local){
    is_sum_based <- all(abs(rowSums(weight_matrix) - 1) < 1e-6)
    if(is_sum_based) {
      emptyset_weights <- rep(1 / num_matrices, num_matrices)
    } else {
      emptyset_weights <- rep(1, num_matrices)
    }
    weight_matrix <- rbind.data.frame(emptyset_weights, weight_matrix)
  }


  # Fonction pour calculer a(e_i)
  calc_a <- function(ei, cluster_i) {
    sum_a <- 0
    for (matrix_index in 1:num_matrices) {
      dist_matrix <- distance_matrices[[matrix_index]]
      if(local){
        weight <- weight_matrix[cluster_i, matrix_index]
      }else{
        weight <- weight_matrix[matrix_index]
      }
      cluster_members <- which(labels == cluster_i)
      sum_distances <- sum(sapply(cluster_members, function(ej) dist_matrix[ei, ej]))
      sum_a <- sum_a + (weight * sum_distances / length(cluster_members))
    }
    return(sum_a)
  }
  
  # Fonction pour calculer b(e_i)
  calc_b <- function(ei, cluster_i) {
    min_b <- Inf
    for (other_cluster in unique(labels)) {
      if (other_cluster == cluster_i) next
      sum_b <- 0
      
      for (matrix_index in 1:num_matrices) {
        dist_matrix <- distance_matrices[[matrix_index]]
        if(local){
          weight <- weight_matrix[cluster_i, matrix_index]
        }else{
          weight <- weight_matrix[matrix_index]
        }
        other_cluster_members <- which(labels == other_cluster)
        sum_distances <- sum(sapply(other_cluster_members, function(ej) dist_matrix[ei, ej]))
        sum_b <- sum_b + (weight * sum_distances / length(other_cluster_members))
      }
      min_b <- min(min_b, sum_b)
    }
    return(min_b)
  }
  
  silhouette_values <- numeric(n)
  for (i in 1:n) {
    cluster_i <- labels[i]
    a_i <- calc_a(i, cluster_i)
    b_i <- calc_b(i, cluster_i)
    silhouette_values[i] <- (b_i - a_i) / max(b_i, a_i)
  }
  
  # Moyenne des silhouettes pour tous les points (ASW)
  asw <- mean(silhouette_values)
  return(list(asw = asw, silhouette_values = silhouette_values))
}




# True libellé
true_lab <- function(y, y_pred) {
  cm <- table(y, y_pred)
  cm_argmax <- apply(cm, 2, which.max)
  y_pred_ <- sapply(y_pred, function(i) cm_argmax[i])
  return(y_pred_)
}


### Data function
transform_table <- function(data, value_name) {
  df <- as.data.frame(data)
  df$id <- 1:nrow(df)
  df %>%
    pivot_longer(cols = -id, names_to = "key", values_to = value_name)
}

### Splite
df_euclu_split <- function(df, type) {
  result <- list()
  
  #meta data
  result$y <- df$y
  result$name <- df$name

  
  if (type == "single") {
    result$x <- df$x
    return(result)
  } else if (type == "multi") {
    df_list <- lapply(1:ncol(df$x), function(i) {
      data.frame(df$x[,i])
    })
    result$x <- df_list
    return(result)
  }
}

### Eucludien
df_euclu_dist <- function(df, type) {
  result <- list()
  
  #meta data
  result$y <- df$y
  result$name <- df$name
  
  #distance
  if (type == "single") {
    result$x <- as.matrix(dist(df$x, method = "euclidean"))
    return(result)
  } else if (type == "multi") {
    result$x <- lapply(df$x, function(col) {
      if(is.numeric(col)) {
        as.matrix(dist(col, method = "euclidean"))
      }
    })
    return(result)
  }
}


# fuzzy to hard
fuzzy_to_hard <- function(clus, vars="clusters"){
  clus <- as.data.frame(clus$clus)
  colnames(clus) <- c("clusters", "membership")
  return(clus[,vars])
}



#one-hote
hamming_dist <- function(a, b) {
  sum(a != b)
}

df_onehot_dist <- function(df, type, dist) {
  
  #meta data
  result <- list()
  result$y <- df$y
  result$name <- df$name
  
  df <- mutate_all(df$x, as.factor)
  
  if (type == 'single' & dist == FALSE) {
    dummies <- dummyVars("~ .", data = df)
    df_onehot <- predict(dummies, newdata = df)
    result$x <- df_onehot
    return(result)
  }
  
  if (type == 'multi' & dist == FALSE) {
    df_list <- lapply(df, function(x) {
      dummies <- dummyVars("~ .", data = data.frame(x))
      df_onehot <- predict(dummies, newdata = data.frame(x))
      return(df_onehot)
    })
    result$x <- df_list
    return(result)
  }
  
  if (type == 'single' & dist == TRUE) {
    dummies <- dummyVars("~ .", data = df)
    df_onehot <- predict(dummies, newdata = df)
    dist_matrix <- as.matrix(dist(df_onehot, method = hamming_dist))
    result$x <- dist_matrix
    return(result)
  }
  
  if (type == 'multi' & dist == TRUE) {
    dist_list <- lapply(df, function(x) {
      dummies <- dummyVars("~ .", data = data.frame(x))
      df_onehot <- predict(dummies, newdata = data.frame(x))
      dist_matrix <- as.matrix(dist(df_onehot, method = hamming_dist))
      return(dist_matrix)
    })
    result$x <- dist_list
    return(result)
  }
}




# multi-view distance
df_multiview_dist <- function(df, type, dist) {
  
  df_list <- df$x
  
  #meta data
  result <- list()
  result$y <- df$y
  result$name <- df$name
  
  if (type == 'single' & dist == FALSE) {
    df_combined <- do.call(cbind, df_list)
    result$x <- df_combined
    return(result)
  }
  
  if (type == 'single' & dist == TRUE) {
    df_combined <- do.call(cbind, df_list)
    dist_matrix <- as.matrix(dist(df_combined))
    result$x <- dist_matrix
    return(result)
  }
  
  if (type == 'multi' & dist == FALSE) {
    result$x <- df_list
    return(result)
  }
  
  if (type == 'multi' & dist == TRUE) {
    dist_list <- lapply(df_list, function(df) {
      dist_matrix <- as.matrix(dist(df))
      return(dist_matrix)
    })
    result$x <- dist_list
    return(result)
  }
}



#OM specification
seqlist <- function(data, list_vars){
  sequences_list <- lapply(list_vars, function(var_name) {
    seqformat(data = data,
              from = "SPELL", to = "STS",
              id = "id", 
              begin = "time", end = "time", 
              status = var_name,
              process = FALSE) %>% 
      seqdef(with.missing = TRUE, right=NA, gaps=NA, left=NA)})
  return(sequences_list)
}


dtw_distance <- function(ts1, ts2) {
  return(dtw(ts1, ts2)$distance)
}

compute_dtw_matrix <- function(df, dimension_col) {
  ids <- unique(df$id)
  n <- length(ids)
  dist_matrix <- matrix(0, n, n)
  rownames(dist_matrix) <- ids
  colnames(dist_matrix) <- ids
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      ts1 <- df %>% filter(id == ids[i]) %>% pull(dimension_col)
      ts2 <- df %>% filter(id == ids[j]) %>% pull(dimension_col)
      dist_matrix[i, j] <- dtw_distance(ts1, ts2)
      dist_matrix[j, i] <- dist_matrix[i, j]
    }
  }
  return(dist_matrix)
}

# time serie
df_timeseries_dist <- function(df, type, method) {
  
  #meta data
  result <- list()
  result$y <- df$y
  result$name <- df$name
  
  df <- df$x
  ts_cols <- grep("^dimension_", names(df), value = TRUE)
  
  if (type == 'multi' && method == 'dtw') {
    dist_matrices <- lapply(ts_cols, function(col) compute_dtw_matrix(df, col))
    names(dist_matrices) <- ts_cols
    
    result$x <-dist_matrices
    return(result)
  }
  
  if (type == 'single' && method == 'dtw') {
    dist_matrices <- lapply(ts_cols, function(col) compute_dtw_matrix(df, col))
    combined_dist_matrix <- Reduce(`+`, dist_matrices)
    
    result$x <-combined_dist_matrix
    return(result)
  }
  
  if (method == 'om') {
    listseq <- seqlist(data=df, list_vars=c("children","married","left" ))
    listdist <- lapply(1:length(listseq), function(i) {
      seqdist(listseq[[i]], method = "OM", sm = "INDELSLOG", with.missing=TRUE)
    })
    
    if(type == 'multi') {
      result$x <- listdist
      return(result)
    }
    
    if (type == 'single') {
      result$x <- Reduce("+", listdist)
      return(result)
    }
  }
}



















### ==>ECM
ECM <- function(df, nclass, beta){
  
  c <- length(unique(df$y))
  yy <- as.integer(df$y)
  
  if(c == nclass){
    ari_ecm <- matrix(NA, nrow = length(beta)*2, ncol = 3 )
    p_ecm <- matrix(NA, nrow = length(beta)*2, ncol = 3 )
    colnames(ari_ecm) <- c('ari','beta','time')
    colnames(p_ecm) <- c('p','beta','time')
    
    ecm_list <- list()
    i <- 0
    for(t in 1:2){
      for(b in beta){
        i <- i+1
        tryCatch({
          ecm <- evclust::ecm(x=df$x, c=c, type='simple',
                              alpha=1,
                              beta=b,
                              delta=10,epsi=1e-3,disp=FALSE)
          y <- as.integer(as.vector(apply(ecm$betp,1,which.max)))
          nclus <- length(unique(y))
        }, error = function(e) {
          assign("nclus", NA, envir = .GlobalEnv)})
        if(!is.na(nclus) & nclus == nclass){
            ari <- extCriteria(yy, y, "Rand")$rand
            prec <- extCriteria(yy, y, "Precision")$precision
          }else{ari <- prec <- ecm <- NA}

        ari_ecm[i,] <- c(ari, b, t)
        p_ecm[i,] <- c(prec, b, t)
        ecm_list[[paste("t", t, "b", b, sep="_")]] <- ecm
        cat('ecm ',i,'\n')
      }
    }
  }else{
    stop("c is different to nclass")
  }
  save(ari_ecm, p_ecm, ecm_list, file = paste0("Result/",df$name,"ecm.Rdata"))
  return(paste('ECM', df$name, 'OK'))
}



##==> ECMDD
ECMDD <- function(df, nclass, beta, delta, gamma, eta){
  
  c <- length(unique(df$y))
  yy <- as.integer(df$y)
  
  if(c == nclass){
    
    delta_values <- delta  
    gamma_values <- gamma
    eta_values <- eta
    
    ari_ecmdd <- matrix(NA, nrow = length(beta)*length(delta_values)*length(gamma_values)*length(eta_values)*2, ncol = 6)
    p_ecmdd <- matrix(NA, nrow = length(beta)*length(delta_values)*length(gamma_values)*length(eta_values)*2, ncol = 6)
    
    colnames(ari_ecmdd) <- c('ari', 'beta', 'delta', 'gamma', 'eta', 'time')
    colnames(p_ecmdd) <- c('p', 'beta', 'delta', 'gamma', 'eta', 'time')
    
    i <- 0
    ecmdd_list <- list()
    for(t in 1:2){
      for(b in beta){
        for(d in delta_values){
          if(d=='Q') dd <- NULL
          for(g in gamma_values){
            for(e in eta_values){
              i <- i + 1
              tryCatch({
                ecmdd <- ecmdd(x=df$x, c=c, type='simple',
                               alpha=1,
                               beta=b,
                               delta=dd, epsi=1e-3, disp=FALSE,
                               gamma=g, eta=e)
  
                y <- as.integer(as.vector(apply(ecmdd$betp,1,which.max)))
                nclus <- length(unique(y))
              }, error = function(e) {
                assign("nclus", NA, envir = .GlobalEnv)})
              if(!is.na(nclus) & nclus == nclass){
                  ari <- extCriteria(yy, y, "Rand")$rand
                  prec <- extCriteria(yy, y, "Precision")$precision
                }else{ari <- prec <- ecmdd <- NA}
              
              ari_ecmdd[i,] <- c(ari, b, d, g, e, t)
              p_ecmdd[i,] <- c(prec, b, d, g, e, t)
              ecmdd_list[[paste("t", t, "b", b, 'd', d, 'g', g, 'e', e, sep="_")]] <- ecmdd
              cat('ecmdd ',i,'\n')
            }
          }
        }
      }
    }
  }else{
    stop("c is different to nclass")
  }
  save(ari_ecmdd, p_ecmdd, ecmdd_list, file = paste0("Result/",df$name,"ecmdd.Rdata"))
  return(paste('ECMDD', df$name, 'OK'))
}


##==> MECMDD_L_S
MECMDD_L_S <- function(df, nclass, beta, delta, gamma, eta, weight='sum', s=NULL){
  
  c <- length(unique(df$y))
  yy <- as.integer(df$y)
  
  if(c == nclass){
    
    delta_values <- delta  
    gamma_values <- gamma
    eta_values <- eta
    
    ari_mecmdd_l_s <- matrix(NA, nrow = length(beta)*length(delta_values)*length(gamma_values)*length(eta_values)*2, ncol = 6)
    p_mecmdd_l_s <- matrix(NA, nrow = length(beta)*length(delta_values)*length(gamma_values)*length(eta_values)*2, ncol = 6)
    
    colnames(ari_mecmdd_l_s) <- c('ari', 'beta', 'delta', 'gamma', 'eta', 'time')
    colnames(p_mecmdd_l_s) <- c('p', 'beta', 'delta', 'gamma', 'eta', 'time')
    
    i <- 0
    mecmdd_l_list <- list()
    for(t in 1:2){
      for(b in beta){
        for(d in delta_values){
          if(d=='Q') dd <- NULL
          for(g in gamma_values){
            for(e in eta_values){
              i <- i + 1
              tryCatch({
                mecmdd_l_s <- mecmdd.rwl(Xlist=df$x, c=c, type='simple',
                                       alpha=1,
                                       beta=b,
                                       delta=dd, epsi=1e-3, disp=FALSE,
                                       gamma=g, eta=e, weight=weight)
                
                y <- as.integer(as.vector(apply(mecmdd_l_s$betp,1,which.max)))
                nclus <- length(unique(y))
              }, error = function(e) {
                assign("nclus", NA, envir = .GlobalEnv)})
              if(!is.na(nclus) & nclus == nclass){
                  ari <- extCriteria(yy, y, "Rand")$rand
                  prec <- extCriteria(yy, y, "Precision")$precision
                }else{ari <- prec <- mecmdd_l_s <- NA}
              
              ari_mecmdd_l_s[i,] <- c(ari, b, d, g, e, t)
              p_mecmdd_l_s[i,] <- c(prec, b, d, g, e, t)
              mecmdd_l_s_list[[paste("t", t, "b", b, 'd', d, 'g', g, 'e', e, sep="_")]] <- mecmdd_l_s
              cat('mecmdd_l_s ',i,'\n')
            }
          }
        }
      }
    }
  }else{
    stop("c is different to nclass")
  }
  save(ari_mecmdd_l_s, p_mecmdd_l_s, mecmdd_l_s_list, file = paste0("Result/",df$name,"mecmdd_l_s.Rdata"))
  return(paste('MECMDD_L_S', df$name, 'OK'))
}


##==> MECMDD_G_S
MECMDD_G_S <- function(df, nclass, beta, delta, gamma, eta, weight='sum', s=NULL){
  
  c <- length(unique(df$y))
  yy <- as.integer(df$y)
  
  if(c == nclass){
    
    delta_values <- delta  
    gamma_values <- gamma
    eta_values <- eta
    
    ari_mecmdd_g_s <- matrix(NA, nrow = length(beta)*length(delta_values)*length(gamma_values)*length(eta_values)*2, ncol = 6)
    p_mecmdd_g_s <- matrix(NA, nrow = length(beta)*length(delta_values)*length(gamma_values)*length(eta_values)*2, ncol = 6)
    
    colnames(ari_mecmdd_g_s) <- c('ari', 'beta', 'delta', 'gamma', 'eta', 'time')
    colnames(p_mecmdd_g_s) <- c('p', 'beta', 'delta', 'gamma', 'eta', 'time')
    
    i <- 0
    mecmdd_g_list <- list()
    for(t in 1:2){
      for(b in beta){
        for(d in delta_values){
          if(d=='Q') dd <- NULL
          for(g in gamma_values){
            for(e in eta_values){
              i <- i + 1
              tryCatch({
                mecmdd_g_s <- mecmdd.rwg(Xlist=df$x, c=c, type='simple',
                                       alpha=1,
                                       beta=b,
                                       delta=dd, epsi=1e-3, disp=FALSE,
                                       gamma=g, eta=e, weight=weight)
                
                y <- as.integer(as.vector(apply(mecmdd_g_s$betp,1,which.max)))
                nclus <- length(unique(y))
              }, error = function(e) {
                assign("nclus", NA, envir = .GlobalEnv)})
              if(!is.na(nclus) & nclus == nclass){
                  ari <- extCriteria(yy, y, "Rand")$rand
                  prec <- extCriteria(yy, y, "Precision")$precision
                }else{ari <- prec <- mecmdd_g_s <- NA}
              
              ari_mecmdd_g_s[i,] <- c(ari, b, d, g, e, t)
              p_mecmdd_g_s[i,] <- c(prec, b, d, g, e, t)
              mecmdd_g_s_list[[paste("t", t, "b", b, 'd', d, 'g', g, 'e', e, sep="_")]] <- mecmdd_g_s
              cat('mecmdd_g_s ',i,'\n')
            }
          }
        }
      }
    }
  }else{
    stop("c is different to nclass")
  }
  save(ari_mecmdd_g_s, p_mecmdd_g_s, mecmdd_g_s_list, file = paste0("Result/",df$name,"mecmdd_g_s.Rdata"))
  return(paste('MECMDD_G_S', df$name, 'OK'))
}





##==> MECMDD_L_P
MECMDD_L_P <- function(df, nclass, beta, delta, gamma, eta, weight='prod', s=NULL){
  
  c <- length(unique(df$y))
  yy <- as.integer(df$y)
  
  if(c == nclass){
    
    delta_values <- delta  
    gamma_values <- gamma
    eta_values <- eta
    
    ari_mecmdd_l_p <- matrix(NA, nrow = length(beta)*length(delta_values)*length(gamma_values)*length(eta_values)*2, ncol = 6)
    p_mecmdd_l_p <- matrix(NA, nrow = length(beta)*length(delta_values)*length(gamma_values)*length(eta_values)*2, ncol = 6)
    
    colnames(ari_mecmdd_l_p) <- c('ari', 'beta', 'delta', 'gamma', 'eta', 'time')
    colnames(p_mecmdd_l_p) <- c('p', 'beta', 'delta', 'gamma', 'eta', 'time')
    
    i <- 0
    mecmdd_l_list <- list()
    for(t in 1:2){
      for(b in beta){
        for(d in delta_values){
          if(d=='Q') dd <- NULL
          for(g in gamma_values){
            for(e in eta_values){
              i <- i + 1
              tryCatch({
                mecmdd_l_p <- mecmdd.rwl(Xlist=df$x, c=c, type='simple',
                                       alpha=1,
                                       beta=b,
                                       delta=dd, epsi=1e-4, disp=FALSE,
                                       gamma=g, eta=e, weight=weight)
                
                y <- as.integer(as.vector(apply(mecmdd_l_p$betp,1,which.max)))
                nclus <- length(unique(y))
              }, error = function(e) {
                assign("nclus", NA, envir = .GlobalEnv)})
              if(!is.na(nclus) & nclus == nclass){
                ari <- extCriteria(yy, y, "Rand")$rand
                prec <- extCriteria(yy, y, "Precision")$precision
              }else{ari <- prec <- mecmdd_l_p <- NA}
              
              ari_mecmdd_l_p[i,] <- c(ari, b, d, g, e, t)
              p_mecmdd_l_p[i,] <- c(prec, b, d, g, e, t)
              mecmdd_l_p_list[[paste("t", t, "b", b, 'd', d, 'g', g, 'e', e, sep="_")]] <- mecmdd_l_p
              cat('mecmdd_l_p ',i ,'\n')
            }
          }
        }
      }
    }
  }else{
    stop("c is different to nclass")
  }
  save(ari_mecmdd_l_p, p_mecmdd_l_p, mecmdd_l_p_list, file = paste0("Result/",df$name,"mecmdd_l_p.Rdata"))
  return(paste('MECMDD_L_P', df$name, 'OK'))
}



##==> MECMDD_G_P
MECMDD_G_P <- function(df, nclass, beta, delta, gamma, eta, weight='prod', s=NULL){
  
  c <- length(unique(df$y))
  yy <- as.integer(df$y)
  
  if(c == nclass){
    
    delta_values <- delta  
    gamma_values <- gamma
    eta_values <- eta
    
    ari_mecmdd_g_p <- matrix(NA, nrow = length(beta)*length(delta_values)*length(gamma_values)*length(eta_values)*2, ncol = 6)
    p_mecmdd_g_p <- matrix(NA, nrow = length(beta)*length(delta_values)*length(gamma_values)*length(eta_values)*2, ncol = 6)
    
    colnames(ari_mecmdd_g_p) <- c('ari', 'beta', 'delta', 'gamma', 'eta', 'time')
    colnames(p_mecmdd_g_p) <- c('p', 'beta', 'delta', 'gamma', 'eta', 'time')
    
    i <- 0
    mecmdd_g_list <- list()
    for(t in 1:2){
      for(b in beta){
        for(d in delta_values){
          if(d=='Q') dd <- NULL
          for(g in gamma_values){
            for(e in eta_values){
              i <- i + 1
              tryCatch({
                mecmdd_g_p <- mecmdd.rwg(Xlist=df$x, c=c, type='simple',
                                       alpha=1,
                                       beta=b,
                                       delta=dd, epsi=1e-4, disp=FALSE,
                                       gamma=g, eta=e, weight=weight)
                
                y <- as.integer(as.vector(apply(mecmdd_g_p$betp,1,which.max)))
                nclus <- length(unique(y))
              }, error = function(e) {
                assign("nclus", NA, envir = .GlobalEnv)})
              if(!is.na(nclus) & nclus == nclass){
                ari <- extCriteria(yy, y, "Rand")$rand
                prec <- extCriteria(yy, y, "Precision")$precision
              }else{ari <- prec <- mecmdd_g_p <- NA}
              
              ari_mecmdd_g_p[i,] <- c(ari, b, d, g, e, t)
              p_mecmdd_g_p[i,] <- c(prec, b, d, g, e, t)
              mecmdd_g_p_list[[paste("t", t, "b", b, 'd', d, 'g', g, 'e', e, sep="_")]] <- mecmdd_g_p
              cat('mecmdd_g_p ',i,'\n')
            }
          }
        }
      }
    }
  }else{
    stop("c is different to nclass")
  }
  save(ari_mecmdd_g_p, p_mecmdd_g_p, mecmdd_g_p_list, file = paste0("Result/",df$name,"mecmdd_g_p.Rdata"))
  return(paste('MECMDD_G_P', df$name, 'OK'))
}




##==> WMVEC
WMVEC <- function(df, nclass, beta, lambda){
  
  c <- length(unique(df$y))
  yy <- as.integer(df$y)
  
  if(c == nclass){
    
    lambda_values <- lambda  
    
    ari_wmvec <- matrix(NA, nrow = length(beta)*length(lambda_values)*2, ncol = 4)
    p_wmvec <- matrix(NA, nrow = length(beta)*length(lambda_values)*2, ncol = 4)
    
    colnames(ari_wmvec) <- c('ari', 'beta', 'lambda', 'time')
    colnames(p_wmvec) <- c('p', 'beta', 'lambda', 'time')
    
    i <- 0
    wmvec_list <- list()
    for(t in 1:2){
      for(b in beta){
        for(l in lambda_values){
          i <- i + 1
          tryCatch({
            wmvec <- wmvec(X=df$x, c=c, alpha = 2, delta = 10, index = 0.5, 
                           n_max = 10, epsilon = 1e-3,  jaccard = 1e10, 
                           J_value = Inf, Normalization = FALSE, 
                           beta = b, lambda = l)
            
            y <- as.integer(as.vector(apply(wmvec$betp,1,which.max)))
            nclus <- length(unique(y))
          }, error = function(e) {
            assign("nclus", NA, envir = .GlobalEnv)})
          if(!is.na(nclus) & nclus == nclass){
              ari <- extCriteria(yy, y, "Rand")$rand
              prec <- extCriteria(yy, y, "Precision")$precision
            }else{ari <- prec <- wmvec <- NA}
          
          ari_wmvec[i,] <- c(ari, b, l, t)
          p_wmvec[i,] <- c(prec, b, l, t)
          wmvec_list[[paste("t", t, "b", b, 'l', l, sep="_")]] <- wmvec
          cat('wmvec ',i,'\n')
        }
      }
    }
  }else{
    stop("c is different to nclass")
  }
  save(ari_wmvec, p_wmvec, wmvec_list, file = paste0("Result/",df$name,"wmvec.Rdata"))
  return(paste('WMVEC', df$name, 'OK'))
}






##==> FCMDD
FCMDD <- function(df, nclass, m){
  
  c <- length(unique(df$y))
  yy <- as.integer(df$y)
  
  if(c == nclass){
    ari_fcmdd <- matrix(NA, nrow = length(m)*2, ncol = 3 )
    p_fcmdd <- matrix(NA, nrow = length(m)*2, ncol = 3 )
    colnames(ari_fcmdd) <- c('ari','m','time')
    colnames(p_fcmdd) <- c('p','m','time')
    
    i <- 0
    fcmdd_list <- list()
    for(t in 1:2){
      for(j in m){
        i <- i+1
        tryCatch({
          fcmdd <- fcmdd(X=df$x, k =c , m = j)
          
          y <- as.integer(as.integer(fuzzy_to_hard(fcmdd, "clusters")))
          nclus <- length(unique(y))
        }, error = function(e) {
          assign("nclus", NA, envir = .GlobalEnv)})
        if(!is.na(nclus) & nclus == nclass){
            ari <- extCriteria(yy, y, "Rand")$rand
            prec <- extCriteria(yy, y, "Precision")$precision
          }else{ari <- prec <- fcmdd <- NA}
        
        ari_fcmdd[i,] <- c(ari, j, t)
        p_fcmdd[i,] <- c(prec, j, t)
        fcmdd_list[[paste("t", t, "m", j, sep="_")]] <- fcmdd
        cat('fcmdd ',i,'\n')
      }
    }
  }else{
    stop("c is different to nclass")
  }
  save(ari_fcmdd, p_fcmdd, fcmdd_list, file = paste0("Result/",df$name,"fcmdd.Rdata"))
  return(paste('FCMDD', df$name, 'OK'))
}





##==> MFCMDD_L
MFCMDD_L <- function(df, nclass, m){
  
  c <- length(unique(df$y))
  yy <- as.integer(df$y)
  
  if(c == nclass){
    ari_mfcmdd_l <- matrix(NA, nrow = length(m)*2, ncol = 3 )
    p_mfcmdd_l <- matrix(NA, nrow = length(m)*2, ncol = 3 )
    colnames(ari_mfcmdd_l) <- c('ari','m','time')
    colnames(p_mfcmdd_l) <- c('p','m','time')
    
    i <- 0
    mfcmdd_l_list <- list()
    for(t in 1:2){
      for(j in m){
        i <- i+1
        tryCatch({
          mfcmdd_l <- mfcmdd.rwl(Xlist=df$x, k=c , m = j,
                                          weight='sum', s=NULL)
          
          y <- as.integer(as.vector(fuzzy_to_hard(mfcmdd_l, "clusters")))
          nclus <- length(unique(y))
        }, error = function(e) {
          assign("nclus", NA, envir = .GlobalEnv)})
        if(!is.na(nclus) & nclus == nclass){
            ari <- extCriteria(yy, y, "Rand")$rand
            prec <- extCriteria(yy, y, "Precision")$precision
          }else{ari <- prec <- mfcmdd_l <- NA}
        
        ari_mfcmdd_l[i,] <- c(ari, j, t)
        p_mfcmdd_l[i,] <- c(prec, j, t)
        mfcmdd_l_list[[paste("t", t, "m", j, sep="_")]] <- mfcmdd_l
        cat('mfcmdd_l ',i,'\n')
      }
    }
  }else{
    stop("c is different to nclass")
  }
  save(ari_mfcmdd_l, p_mfcmdd_l, mfcmdd_l_list, file = paste0("Result/",df$name,"mfcmdd_l.Rdata"))
  return(paste('MFCMDD_L', df$name, 'OK'))
}



##==> MFCMDD_G
MFCMDD_G <- function(df, nclass, m){
  
  c <- length(unique(df$y))
  yy <- as.integer(df$y)
  
  if(c == nclass){
    ari_mfcmdd_g <- matrix(NA, nrow = length(m)*2, ncol = 3 )
    p_mfcmdd_g <- matrix(NA, nrow = length(m)*2, ncol = 3 )
    colnames(ari_mfcmdd_g) <- c('ari','m','time')
    colnames(p_mfcmdd_g) <- c('p','m','time')
    
    i <- 0
    mfcmdd_g_list <- list()
    for(t in 1:2){
      for(j in m){
        i <- i+1
        tryCatch({
          mfcmdd_g <- mfcmdd.rwg(Xlist=df$x, k=c , m = j,
                                    weight='sum', s=NULL)
          
          y <- as.integer(as.vector(fuzzy_to_hard(mfcmdd_g, "clusters")))
          nclus <- length(unique(y))
        }, error = function(e) {
          assign("nclus", NA, envir = .GlobalEnv)})
        if(!is.na(nclus) & nclus == nclass){
            ari <- extCriteria(yy, y, "Rand")$rand
            prec <- extCriteria(yy, y, "Precision")$precision
          }else{ari <- prec <- mfcmdd_g <- NA}
        
        ari_mfcmdd_g[i,] <- c(ari, j, t)
        p_mfcmdd_g[i,] <- c(prec, j, t)
        mfcmdd_g_list[[paste("t", t, "m", j, sep="_")]] <- mfcmdd_g
        cat('mfcmdd_g ',i,'\n')
      }
    }
  }else{
    stop("c is different to nclass")
  }
  save(ari_mfcmdd_g, p_mfcmdd_g, mfcmdd_g_list, file = paste0("Result/",df$name,"mfcmdd_g.Rdata"))
  return(paste('MFCMDD_G', df$name, 'OK'))
}




##==> MRDCA_L
MRDCA_L<- function(df, nclass){
  
  c <- length(unique(df$y))
  yy <- as.integer(df$y)
  
  if(c == nclass){
    ari_mrdca_l <- matrix(NA, nrow = 2, ncol = 2 )
    p_mrdca_l <- matrix(NA, nrow = 2, ncol = 2 )
    colnames(ari_mrdca_l) <- c('ari','time')
    colnames(p_mrdca_l) <- c('p','time')
    
    mrdca_l_list <- list()
    i <- 0
    for(t in 1:2){
      i <- i+1
      tryCatch({
        mrdca_l <- mrdca.rwl(k=c, q=1, data=df$x)
          
        y <- as.integer(as.vector(apply(mrdca_l$clus,1,which.max)))
        nclus <- length(unique(y))
      }, error = function(e) {
        assign("nclus", NA, envir = .GlobalEnv)})
      if(!is.na(nclus) & nclus == nclass){
          ari <- extCriteria(yy, y, "Rand")$rand
          prec <- extCriteria(yy, y, "Precision")$precision
        }else{ari <- prec <- mrdca_l <- NA}
      
      ari_mrdca_l[i,] <- c(ari, t)
      p_mrdca_l[i,] <- c(prec, t)
      mrdca_l_list[[paste("t", t, sep="_")]] <- mrdca_l
      cat('mrdca_l ',i,'\n')
    }
  }else{
    stop("c is different to nclass")
  }
  save(ari_mrdca_l, p_mrdca_l, mrdca_l_list, file = paste0("Result/",df$name,"mrdca_l.Rdata"))
  return(paste('MRDCA_L', df$name, 'OK'))
}



##==> MRDCA_G
MRDCA_G<- function(df, nclass){
  
  c <- length(unique(df$y))
  yy <- as.integer(df$y)
  
  if(c == nclass){
    ari_mrdca_g <- matrix(NA, nrow = 2, ncol = 2 )
    p_mrdca_g <- matrix(NA, nrow = 2, ncol = 2 )
    colnames(ari_mrdca_g) <- c('ari','time')
    colnames(p_mrdca_g) <- c('p','time')
    
    mrdca_g_list <- list()
    i <- 0
    for(t in 1:2){
      i <- i+1
      tryCatch({
        mrdca_g <- mrdca.rwg(k=c, q=1, data=df$x)
        
        y <- as.integer(as.vector(apply(mrdca_g$clus,1,which.max)))
        nclus <- length(unique(y))
      }, error = function(e) {
        assign("nclus", NA, envir = .GlobalEnv)})
      if(!is.na(nclus) & nclus == nclass){
        ari <- extCriteria(yy, y, "Rand")$rand
        prec <- extCriteria(yy, y, "Precision")$precision
      }else{ari <- prec <- mrdca_g <- NA}
      
      ari_mrdca_g[i,] <- c(ari, t)
      p_mrdca_g[i,] <- c(prec, t)
      mrdca_g_list[[paste("t", t, sep="_")]] <- mrdca_g
      cat('mrdca_g ',i)
    }
  }else{
    stop("c is different to nclass")
  }
  save(ari_mrdca_g, p_mrdca_g, mrdca_g_list, file = paste0("Result/",df$name,"mrdca_g.Rdata"))
  return(paste('MRDCA_G', df$name, 'OK'))
}