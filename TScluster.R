library('TSdist')  # for TSDistances
library('tseries') # for ar
library('MASS')    # for isoMDS
library('Kmedians')# for Kmedians

# create data

phi <- c(0.25, 0.5)
var <- 1
n <- 100
for(i in 1:num_series){
  ar <- arima.sim(n = n, list(ar = phi), sd = sqrt(var))
  data <- rbind(data, ar)
}

# create dissimilarity matrix

distance_matrix <- matrix(NA, nrow = nrow(data), ncol = nrow(data))
for(i in 1:nrow(data)){
  for(j in 1:nrow(data)){
    distance_matrix[i, j] <- TSDistances(data[i,], data[j,], distance = "ar.lpc.ceps")
  }
}

# find point through Kruskal

data_new <- isoMDS(distance_matrix)$points

# use Kruskal to cluster

kmeans(data_new, 3)

Kmedians(data_new, 3)

# create cepstrum

cepstra <- function(x, n){
  c <- c()
  for(i in 1:n){
    if (i == 1) {
      c[i] <- x[i]
    }else if(1 < i && i <= length(x)){
      param <- 0
      for(j in 1:(i - 1)){
        param <- (1 - (j/i))*x[j]*c[i - j] + param
      }
      c[i] <- -x[i] - param
    }else if(i > length(x)){
      param2 <- 0
      for(j in 1:length(x)){
        param2 <- (1 - j/i)*x[j]*c[i - j] + param2
      }
      c[i] <- -param2
    }
  }
  return(c)
}
cep <- matrix(NA, nrow = nrow(data), ncol = 10)
for(i in 1:nrow(data)){
  cep[i, ] <- cepstra(ar(data[i, ])$ar, 10)
}

# use cepstra to do cluster

kmeans(cep, 3)

Kmedians(cep, 3)
