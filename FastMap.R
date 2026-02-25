FastMap <- function(D, k){
  n <- nrow(D)
  X <- matrix(0, nrow = n, ncol = k)
  new_D <- D
  for(l in 1:k){
    D_2 <- new_D^2
    b <- sample(1:n, 1)
    for(i in 1:5){
      a <- which.max(new_D[b, ])
      b <- which.max(new_D[a, ])
    }
    if(a == b){
      break
    }
    for(i in 1:n){
      X[i, l] <- (D_2[a, i] + D_2[a, b] - D_2[b, i])/(2*new_D[a, b])  
    }
    for(i in 1:n){
      for(j in 1:i){
        new_D[i, j] <- sqrt(max(D_2[i, j] - (X[i, l] - X[j, l])^2, 0))
        new_D[j, i] <- new_D[i, j]
      }
    }
    if(max(new_D) < 1e-6){
      break
    }
  }
  return(X)
}

Sim <- function(Ground, Cluster){
  G_groups <- unique(Ground)
  A_groups <- unique(Cluster)
  
  k <- length(G_groups)
  sim_sum <- 0
  
  for (g in G_groups) {
    G_i_idx <- which(Ground == g)
    
    max_sim <- 0
    for (a in A_groups) {
      A_j_idx <- which(Cluster == a)
      
      intersect_size <- length(intersect(G_i_idx, A_j_idx))
      sim_ij <- (2 * intersect_size) / (length(G_i_idx) + length(A_j_idx))
      
      if (sim_ij > max_sim) {
        max_sim <- sim_ij
      }
    }
    
    sim_sum <- sim_sum + max_sim
  }
  
  return(sim_sum / k)
}

# 模擬資料

set.seed(1)
data <- matrix(rnorm(20), ncol = 2)
D <- as.matrix(dist(data, method = "euclidean"))
X <- FastMap(D, k = 2)