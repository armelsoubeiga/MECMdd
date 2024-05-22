
# %%%
update_prototype <- function(P, weight, data) {

  new_prototypes <- numeric(length(P))
  for (k in 1:length(P)) {
    sum_distances <- sapply(1:nrow(data[[1]]), function(i) {
      sum(sapply(1:length(data), function(j) {
        weight[k,j] * data[[j]][i, P[[k]]]
      }))
    })
    new_prototypes[k] <- which.min(sum_distances)
  }
  return(new_prototypes)
}


update_prototype_global <- function(P, weight, data) {
  
  new_prototypes <- numeric(length(P))
  for (k in 1:length(P)) {
    sum_distances <- sapply(1:nrow(data[[1]]), function(i) {
      sum(sapply(1:length(data), function(j) {
        weight[j] * data[[j]][i, P[[k]]]
      }))
    })
    new_prototypes[k] <- which.min(sum_distances)
  }
  return(new_prototypes)
}




# %%%
update_weight_local <- function(P, prototype, data) {
  
  p <- length(data)
  new_weights <- matrix(0, nrow = length(P), ncol = length(data))
  for (k in 1:length(P)) {
    sum_distances <- sapply(1:length(data), function(j) {
      sum(data[[j]][P[[k]], prototype[k]])
    })
    prod_sum_distances <- prod(sum_distances)^(1/p)
    new_weights[k,] <- prod_sum_distances / sum_distances
  }

  return(new_weights)
}


update_weight_global <- function(P, prototype, data) {

  p <- length(data)
  new_weights <- numeric(length(data))
  
  for (j in 1:length(data)) {
    sum_distances <- sapply(1:length(P), function(k) {
      sum(data[[j]][P[[k]], prototype[k]])
    })
    prod_sum_distances <- prod(sapply(1:length(data), function(h) {
      sum(sapply(1:length(P), function(k) {
        sum(data[[h]][P[[k]], prototype[k]])
      }))
    }))^(1/p)
    new_weights[j] <- prod_sum_distances / sum(sum_distances)
  }
  return(new_weights)
}



# %%
update_partition <- function(P, prototype, weight, data) {

  J <- 0
  for (k in 1:length(P)) {
    sum_distances <- sum(sapply(1:length(data), function(j) {
      weight[k,j] * sum(data[[j]][P[[k]], prototype[k]])
    }))

    J <- J + sum_distances
  }
  return(J)
}


update_partition_global <- function(P, prototype, weight, data) {
  
  J <- 0
  for (k in 1:length(P)) {
    sum_distances <- sum(sapply(1:length(data), function(j) {
      weight[j] * sum(data[[j]][P[[k]], prototype[k]])
    }))
    
    J <- J + sum_distances
  }
  return(J)
}




# %%
cluster_index <- function(n,P){
  clus <- matrix(0, nrow = n, ncol = length(P))
  for (k in 1:length(P)) {
    for (i in P[[k]]) {
      clus[i, k] <- 1
    }
  }
  return(clus)
}
