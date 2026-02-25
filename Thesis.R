install.packages(c("sf", "dplyr", "ggplot2", "rnaturalearth", "rnaturalearthdata"))
library('sf')
library('dplyr')
library('ggplot2')
library('rnaturalearth')
library('rnaturalearthdata')
library('TSdist')
library('cluster')
library('aTSA')
library('forecast')
library('ClustGeo')
library('xtable')
library('TTR')
library('zoo')
source("D:/Code/R/FastMap.R")
source("D:/Code/R/TScluster.R")

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
      c[i] <- x[i] + param
    }else if(i > length(x)){
      param2 <- 0
      for(j in 1:length(x)){
        param2 <- (1 - j/i)*x[j]*c[i - j] + param2
      }
      c[i] <- param2
    }
  }
  return(c)
}

folder_path <- "D:/Code/IncomeData"  

file_list <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)

data_list <- lapply(file_list, function(file) {
  read.table(file, header = FALSE)  
})

combined_data <- do.call(cbind, data_list)

colnames(combined_data) <- tools::file_path_sans_ext(basename(file_list))

for(i in 1:ncol(combined_data)){
  print(ndiffs(combined_data[, i]))
}

combined_data_new <- as.data.frame(
  apply(combined_data, 2, function(x) rollmean(x, k = 2, align = 'center'))
)

log_transformed_data <- as.data.frame(lapply(combined_data_new, log))

ar <- c()
for(i in 1:ncol(log_transformed_data)){
  ar[i] <- arima(log_transformed_data[, i], order = c(1,1,0))$coef
}

cep <- matrix(0, nrow = length(ar), ncol = 10)
for(i in 1:length(ar)){
  cep[i, ] <- cepstra(ar[i], 10)
}

LPC <- dist(cep)
LPC <- (LPC - min(LPC))/(max(LPC)-min(LPC))

cluster_kmedoids_lpc_or <- pam(LPC, 2, diss= T)$clustering
mean(silhouette(cluster_kmedoids_lpc_or, LPC)[, 3])
Sim(ground_truth, cluster_kmedoids_lpc_or)

series_names <- colnames(combined_data)

n <- ncol(combined_data)

dist_matrix <- matrix(0, nrow = n, ncol = n)
rownames(dist_matrix) <- series_names
colnames(dist_matrix) <- series_names

for (i in 1:n) {
  for (j in i:n) {  
    dist_ij <- TSDistances(log_transformed_data[[i]], log_transformed_data[[j]], distance = "cdm")
    dist_matrix[i, j] <- dist_ij
    dist_matrix[j, i] <- dist_ij
  }
}

CDM <- dist_matrix
CDM <- (CDM - min(CDM))/(max(CDM)-min(CDM))
diag(CDM) <- 0

coords_latlon <- as.matrix(ordered_coords[, c("lon", "lat")])
euclidean_dist_matrix <- as.matrix(dist(coords_latlon, method = "euclidean"))

rownames(euclidean_dist_matrix) <- states_csv_sorted$state
colnames(euclidean_dist_matrix) <- states_csv_sorted$state

GEO <- (euclidean_dist_matrix - min(euclidean_dist_matrix))/(max(euclidean_dist_matrix)-min(euclidean_dist_matrix))

datap_cdm <- FastMap(CDM, 24)
datap_lpc <- FastMap(as.matrix(LPC), 2)
coord <- FastMap(GEO, 2)
LPC <- as.dist(LPC)
CDM <- as.dist(CDM)

alpha <- seq(0, 1, by = 0.1)

weight <- 0.04*0.04/2

Q0_t_cdm <- c()
Q0_t_lpc <- c()
Q1_t <- c()
Q0_t_cdm <- weight*sum(dist(datap_cdm)**2)
Q0_t_lpc <- weight*sum(dist(datap_lpc)**2)
Q1_t <- weight*sum(dist(coord)**2)

# test best K

SSE_test <- c()
for(K in 1:10){
  clus <- kmeans(datap_lpc, K)
  SSE_test[K] <- clus$tot.withinss
}
plot(1:10, SSE_test, xlab = 'K', ylab = 'SSE')

alpha_test <- seq(0, 1, 0.1)
cluster_kmean_lpc_test <- matrix(0, nrow = n, ncol = length(alpha_test))
SSE_lpc_test <- c()
for(i in 1:length(alpha_test)){
  data_lpc <- cbind(sqrt(weight*2*(1-alpha_test[i])/0.08)*datap_lpc, sqrt(weight*2*alpha_test[i]/0.08)*coord)
  clustering_kmean_test <- kmeans(data_lpc, 3)
  cluster_kmean_lpc_test[,i] <- clustering_kmean_test$cluster
  SSE_lpc_test[i]<- clustering_kmean_test$tot.withinss
}

Q0_i_lpc_test <- matrix(0, nrow = 3, ncol = length(alpha_test))
for(i in 1:length(alpha_test)){
  for(j in 1:3){
    if(sum(cluster_kmean_lpc_test[,i] == j) == 1){
      Q0_i_lpc_test [j, i] <- 0
    }else{
      Q0_i_lpc_test [j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmean_lpc_test[,i] == j))*sum(dist(datap_lpc[which(cluster_kmean_lpc_test[,i] == j), ])**2)
    }
  }
}
Q1_i_lpc_test  <- matrix(0, nrow = 3, ncol = length(alpha_test))
for(i in 1:length(alpha_test)){
  for(j in 1:3){
    if(sum(cluster_kmean_lpc_test[,i] == j) == 1){
      Q1_i_lpc_test [j, i] <- 0
    }else{
      Q1_i_lpc_test [j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmean_lpc_test[,i] == j))*sum(dist(coord[which(cluster_kmean_lpc_test[,i] == j), ])**2)
    }
  }
}
1-(colSums(Q0_i_lpc_test)/Q0_t_lpc)
1-(colSums(Q1_i_lpc_test)/Q1_t)

## K-means
SSE_test <- c()
for(K in 1:10){
  clus <- kmeans(datap_lpc, K)
  SSE_test[K] <- clus$tot.withinss
}
plot(1:10, SSE_test, xlab = 'K', ylab = 'SSE')

k_vals <- c(1:10)
points <- cbind(k_vals, SSE_test)

A <- points[1, ]
B <- points[nrow(points), ]

point_to_line_dist <- function(P, A, B) {
  num <- abs((B[2] - A[2]) * P[1] - (B[1] - A[1]) * P[2] + B[1]*A[2] - B[2]*A[1])
  denom <- sqrt((B[2] - A[2])^2 + (B[1] - A[1])^2)
  return(num / denom)
}

dists <- sapply(2:(nrow(points)-1), function(i) point_to_line_dist(points[i, ], A, B))

elbow_index <- which.max(dists) + 1
elbow_k <- k_vals[elbow_index]

cat("最佳分群數 K =", elbow_k, "\n")

## K-medoids

SSE_test <- c()
for(K in 1:10){
  clus <- pam(LPC, K, diss = T)
  SSE_test[K] <- clus$objective[2]
}

k_vals <- c(1:10)
plot(k_vals, SSE_test, xlab = 'K', ylab = 'S')

points <- cbind(k_vals, SSE_test)

A <- points[1, ]
B <- points[nrow(points), ]

point_to_line_dist <- function(P, A, B) {
  num <- abs((B[2] - A[2]) * P[1] - (B[1] - A[1]) * P[2] + B[1]*A[2] - B[2]*A[1])
  denom <- sqrt((B[2] - A[2])^2 + (B[1] - A[1])^2)
  return(num / denom)
}

dists <- sapply(2:(nrow(points)-1), function(i) point_to_line_dist(points[i, ], A, B))

elbow_index <- which.max(dists) + 1
elbow_k <- k_vals[elbow_index]

cat("最佳分群數 K =", elbow_k, "\n")

## Hierarchical
tree <- hclustgeo(as.dist(LPC))
plot(tree,hang=-1, label=FALSE, xlab="", sub="", main="")


us_states_selected <- us_states_clustered %>%
  filter(!is.na(state))

us_states_selected <- us_states_selected %>%
  mutate(centroid = st_centroid(geometry))

coords <- st_coordinates(us_states_selected$centroid)
us_states_selected <- us_states_selected %>%
  mutate(
    lon = coords[, 1],
    lat = coords[, 2]
  )
state_coords <- us_states_selected %>%
  dplyr::select(state, lon, lat)
ordered_coords <- us_states_selected %>%
  filter(state %in% series_names) %>%
  arrange(match(state, series_names)) %>%
  dplyr::select(state, lon, lat) %>%
  st_drop_geometry()


ts_data <- ts(combined_data[, 25], start = 1929, end = 1999, frequency = 1)

plot(ts_data, main = colnames(combined_data)[25], xlab = "Year", ylab = "Value",
     col = "black", type = "l")

length(colnames(combined_data))

# Cluster

## K-means

### CDM


cluster_kmean_cdm <- matrix(0, nrow = n, ncol = length(alpha))
SSE_cdm <- c()
for(i in 1:length(alpha)){
  data_cdm <- cbind(sqrt(weight*2*(1-alpha[i])/0.08)*datap_cdm, sqrt(weight*2*alpha[i]/0.08)*coord)
  clustering_kmean <- kmeans(data_cdm, 2)
  cluster_kmean_cdm[,i] <- clustering_kmean$cluster
  SSE_cdm[i]<- clustering_kmean$tot.withinss
}
Q0_i_cdm <- matrix(0, nrow = 2, ncol = length(alpha))
for(i in 1:length(alpha)){
  for(j in 1:2){
    Q0_i_cdm[j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmean_cdm[,i] == j))*sum(dist(datap_cdm[which(cluster_kmean_cdm[,i] == j), ])**2)
  }
}
Q1_i_cdm <- matrix(0, nrow = 2, ncol = length(alpha))
for(i in 1:length(alpha)){
  for(j in 1:2){
    Q1_i_cdm[j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmean_cdm[,i] == j))*sum(dist(coord[which(cluster_kmean_cdm[,i] == j), ])**2)
  }
}

Q0 <- 1-(colSums(Q0_i_cdm)/Q0_t_cdm)
Q1 <- 1-(colSums(Q1_i_cdm)/Q1_t)
1-(colSums(Q0_i_cdm)/Q0_t_cdm)+1-(colSums(Q1_i_cdm)/Q1_t)
SSE_cdm

mat <- matrix(c(Q0, Q1, Q0+Q1, SSE_cdm), nrow = length(Q0), byrow = F)
latex_code <- xtable(mat, digits = c(0, rep(3, ncol(mat))))
print(latex_code, type = "latex")

cluster_kmeans_cdm_or <- kmeans(datap_cdm, 2)$cluster
mean(silhouette(cluster_kmeans_cdm_or, CDM)[, 3])
mean(silhouette(cluster_kmean_cdm[,10], sqrt(0.02*(1-0.9)*CDM**2 + 0.02*0.9*GEO**2))[, 3])
Sim(ground_truth, cluster_kmeans_cdm_or)
Sim(ground_truth, cluster_kmean_cdm[, 10])
### LPC


cluster_kmean_lpc <- matrix(0, nrow = n, ncol = length(alpha))
SSE_lpc <- c()
for(i in 1:length(alpha)){
  data_lpc <- cbind(sqrt(weight*2*(1-alpha[i])/0.08)*datap_lpc, sqrt(weight*2*alpha[i]/0.08)*coord)
  clustering_kmean <- kmeans(data_lpc, 2)
  cluster_kmean_lpc[,i] <- clustering_kmean$cluster
  SSE_lpc[i]<- clustering_kmean$tot.withinss
}

Q0_i_lpc <- matrix(0, nrow = 2, ncol = length(alpha))
for(i in 1:length(alpha)){
  for(j in 1:2){
    Q0_i_lpc[j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmean_lpc[,i] == j))*sum(dist(datap_lpc[which(cluster_kmean_lpc[,i] == j), ])**2)
  }
}
Q1_i_lpc <- matrix(0, nrow = 2, ncol = length(alpha))
for(i in 1:length(alpha)){
  for(j in 1:2){
    Q1_i_lpc[j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmean_lpc[,i] == j))*sum(dist(coord[which(cluster_kmean_lpc[,i] == j), ])**2)
  }
}
Q0 <- 1-(colSums(Q0_i_lpc)/Q0_t_lpc)
Q1 <- 1-(colSums(Q1_i_lpc)/Q1_t)
1-(colSums(Q0_i_lpc)/Q0_t_lpc) + 1-(colSums(Q1_i_lpc)/Q1_t)
SSE_lpc[1:11]

mat <- matrix(c(Q0, Q1, Q0+Q1, SSE_lpc), nrow = length(Q0), byrow = F)
latex_code <- xtable(mat, digits = c(0, rep(4, ncol(mat))))
print(latex_code, type = "latex")


cluster_kmeans_lpc_or <- kmeans(datap_lpc, 2)$cluster
mean(silhouette(cluster_kmeans_lpc_or, LPC)[, 3])
mean(silhouette(cluster_kmean_lpc[,6], sqrt(0.02*(1-0.5)*as.matrix(LPC)**2 + 0.02*0.5*GEO**2))[, 3])
Sim(ground_truth, cluster_kmeans_lpc_or)
Sim(ground_truth, cluster_kmean_lpc[, 6])

### K = 3

cluster_kmean_lpc <- matrix(0, nrow = n, ncol = length(alpha))
SSE_lpc <- c()
for(i in 1:length(alpha)){
  data_lpc <- cbind(sqrt(weight*2*(1-alpha[i])/0.08)*datap_lpc, sqrt(weight*2*alpha[i]/0.08)*coord)
  clustering_kmean <- kmeans(data_lpc, 3)
  cluster_kmean_lpc[,i] <- clustering_kmean$cluster
  SSE_lpc[i]<- clustering_kmean$tot.withinss
}

Q0_i_lpc <- matrix(0, nrow = 3, ncol = length(alpha))
for(i in 1:length(alpha)){
  for(j in 1:3){
    Q0_i_lpc[j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmean_lpc[,i] == j))*sum(dist(datap_lpc[which(cluster_kmean_lpc[,i] == j), ])**2)
  }
}
Q1_i_lpc <- matrix(0, nrow = 3, ncol = length(alpha))
for(i in 1:length(alpha)){
  for(j in 1:3){
    Q1_i_lpc[j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmean_lpc[,i] == j))*sum(dist(coord[which(cluster_kmean_lpc[,i] == j), ])**2)
  }
}
Q0 <- 1-(colSums(Q0_i_lpc)/Q0_t_lpc)
Q1 <- 1-(colSums(Q1_i_lpc)/Q1_t)

mat <- matrix(c(Q0, Q1, Q0+Q1, SSE_lpc), nrow = length(Q0), byrow = F)
latex_code <- xtable(mat, digits = c(0, rep(3, ncol(mat))))
print(latex_code, type = "latex")

mean(silhouette(cluster_kmean_lpc[,4], sqrt(0.02*(1-0.3)*as.matrix(LPC)**2 + 0.02*0.3*GEO**2))[, 3])
Sim(ground_truth_new, cluster_kmean_lpc[, 4])


## K-medoids

### CDM
cluster_kmedoids_cdm <- matrix(0, nrow = n, ncol = length(alpha))
dis_cdm <- c()
for(i in 1:length(alpha)){
  D_CDM <- sqrt(0.02*(1-alpha[i])*CDM**2 + 0.02*alpha[i]*GEO**2)
  k_medoids <- pam(D_CDM, 2, diss= T)
  cluster_kmedoids_cdm[, i] <- k_medoids$clustering
  dis_cdm[i] <- k_medoids$objective[2]
}
Q0_i_cdm_medoid <- matrix(0, nrow = 2, ncol = length(alpha))
for(i in 1:length(alpha)){
  for(j in 1:2){
    if(sum(cluster_kmedoids_cdm[,i] == j) == 1){
      Q0_i_cdm_medoid[j, i] <- 0
    }else{
      Q0_i_cdm_medoid[j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmedoids_cdm[,i] == j))*sum(dist(datap_cdm[which(cluster_kmedoids_cdm[,i] == j), ])**2)
    }
  }
}
Q1_i_cdm_medoid <- matrix(0, nrow = 2, ncol = length(alpha))
for(i in 1:length(alpha)){
  for(j in 1:2){
    if(sum(cluster_kmedoids_cdm[,i] == j) == 1){
      Q1_i_cdm_medoid[j, i] <- 0
    }else{
      Q1_i_cdm_medoid[j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmedoids_cdm[,i] == j))*sum(dist(coord[which(cluster_kmedoids_cdm[,i] == j), ])**2)
    }
  }
}

Q0 <- 1-(colSums(Q0_i_cdm_medoid)/Q0_t_cdm)
Q1 <- 1-(colSums(Q1_i_cdm_medoid)/Q1_t)
mat <- matrix(c(Q0, Q1, Q0+Q1, dis_cdm), nrow = length(Q0), byrow = F)
latex_code <- xtable(mat, digits = c(0, rep(3, ncol(mat))))
print(latex_code, type = "latex")

cluster_kmedoids_cdm_or <- pam(CDM, 2, diss= T)$clustering
mean(silhouette(cluster_kmedoids_cdm_or, CDM)[, 3])
mean(silhouette(cluster_kmedoids_cdm[,5], sqrt(0.02*(1-0.4)*CDM**2 + 0.02*0.4*GEO**2))[, 3])
Sim(ground_truth, cluster_kmedoids_cdm_or)
Sim(ground_truth, cluster_kmedoids_cdm[, 4])

### LPC
cluster_kmedoids_lpc <- matrix(0, nrow = n, ncol = length(alpha))
dis_lpc <- c()
for(i in 1:length(alpha)){
  D_LPC <- sqrt(0.02*(1-alpha[i])*as.matrix(LPC)**2 + 0.02*alpha[i]*GEO**2)
  k_medoids <- pam(D_LPC, 2, diss= T)
  cluster_kmedoids_lpc[, i] <- k_medoids$clustering
  dis_lpc[i] <- k_medoids$objective[2]
}
Q0_i_lpc_medoid <- matrix(0, nrow = 2, ncol = length(alpha))
for(i in 1:length(alpha)){
  for(j in 1:2){
    if(sum(cluster_kmedoids_lpc[,i] == j) == 1){
      Q0_i_lpc_medoid[j, i] <- 0
    }else{
      Q0_i_lpc_medoid[j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmedoids_lpc[,i] == j))*sum(dist(datap_lpc[which(cluster_kmedoids_lpc[,i] == j), ])**2)
    }
  }
}
Q1_i_lpc_medoid <- matrix(0, nrow = 2, ncol = length(alpha))
for(i in 1:length(alpha)){
  for(j in 1:2){
    if(sum(cluster_kmedoids_lpc[,i] == j) == 1){
      Q1_i_lpc_medoid[j, i] <- 0
    }else{
      Q1_i_lpc_medoid[j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmedoids_lpc[,i] == j))*sum(dist(coord[which(cluster_kmedoids_lpc[,i] == j), ])**2)
    }
  }
}

Q0 <- 1-(colSums(Q0_i_lpc_medoid)/Q0_t_lpc)
Q1 <- 1-(colSums(Q1_i_lpc_medoid)/Q1_t)
mat <- matrix(c(Q0, Q1, Q0+Q1, dis_lpc), nrow = length(Q0), byrow = F)
latex_code <- xtable(mat, digits = c(0, rep(5, ncol(mat))))
print(latex_code, type = "latex")

cluster_kmedoids_lpc_or <- pam(LPC, 2, diss= T)$clustering
mean(silhouette(cluster_kmedoids_lpc_or, LPC)[, 3])
mean(silhouette(cluster_kmedoids_lpc[,5], sqrt(0.02*(1-0.4)*as.matrix(LPC)**2 + 0.02*0.4*GEO**2))[, 3])
Sim(ground_truth, cluster_kmedoids_lpc_or)
Sim(ground_truth, cluster_kmedoids_lpc[, 5])
?pam

### K = 3
cluster_kmedoids_lpc <- matrix(0, nrow = n, ncol = length(alpha))
dis_lpc <- c()
for(i in 1:length(alpha)){
  D_LPC <- sqrt(0.02*(1-alpha[i])*as.matrix(LPC)**2 + 0.02*alpha[i]*GEO**2)
  k_medoids <- pam(D_LPC, 3, diss= T)
  cluster_kmedoids_lpc[, i] <- k_medoids$clustering
  dis_lpc[i] <- k_medoids$objective[2]
}
Q0_i_lpc_medoid <- matrix(0, nrow = 3, ncol = length(alpha))
for(i in 1:length(alpha)){
  for(j in 1:3){
    if(sum(cluster_kmedoids_lpc[,i] == j) == 1){
      Q0_i_lpc_medoid[j, i] <- 0
    }else{
      Q0_i_lpc_medoid[j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmedoids_lpc[,i] == j))*sum(dist(datap_lpc[which(cluster_kmedoids_lpc[,i] == j), ])**2)
    }
  }
}
Q1_i_lpc_medoid <- matrix(0, nrow = 3, ncol = length(alpha))
for(i in 1:length(alpha)){
  for(j in 1:3){
    if(sum(cluster_kmedoids_lpc[,i] == j) == 1){
      Q1_i_lpc_medoid[j, i] <- 0
    }else{
      Q1_i_lpc_medoid[j, i] <- 0.04*0.04/(2*0.04*sum(cluster_kmedoids_lpc[,i] == j))*sum(dist(coord[which(cluster_kmedoids_lpc[,i] == j), ])**2)
    }
  }
}

Q0 <- 1-(colSums(Q0_i_lpc_medoid)/Q0_t_lpc)
Q1 <- 1-(colSums(Q1_i_lpc_medoid)/Q1_t)
mat <- matrix(c(Q0, Q1, Q0+Q1, dis_lpc), nrow = length(Q0), byrow = F)
latex_code <- xtable(mat, digits = c(0, rep(3, ncol(mat))))
print(latex_code, type = "latex")

mean(silhouette(cluster_kmedoids_lpc[,9], sqrt(0.02*(1-0.8)*as.matrix(LPC)**2 + 0.02*0.8*GEO**2))[, 3])
Sim(ground_truth_new, cluster_kmedoids_lpc[, 9])

## Hierarchical
### CDM

choicealpha <- choicealpha(as.dist(CDM), as.dist(GEO), range.alpha = alpha, K = 2, scale = F)
sil_cdm <- c()
I_cdm <- c()
for(i in 1:length(alpha)){
  D_CDM <- sqrt(0.02*(1-alpha[i])*CDM**2 + 0.02*alpha[i]*GEO**2)
  cdm_tree <- hclustgeo(as.dist(CDM), as.dist(GEO), alpha = alpha[i])
  cdm_hie <- cutree(cdm_tree, 2)
  I_cdm[i] <- sum(cdm_tree$height[1:23]) 
  sil_cdm[i] <- mean(silhouette(cdm_hie, D_CDM)[, 3])
}

Q0 <- choicealpha$Q[,1]
Q1 <- choicealpha$Q[,2]
mat <- matrix(c(Q0, Q1, Q0+Q1, sil_cdm), nrow = length(Q0), byrow = F)
latex_code <- xtable(mat, digits = c(0, rep(3, ncol(mat))))
print(latex_code, type = "latex")

tree <- hclustgeo(as.dist(CDM))
cdm_hie_org <- cutree(tree, 2)
cdm_tree <- hclustgeo(as.dist(CDM), as.dist(GEO), alpha = 0.1)
cdm_hie <- cutree(cdm_tree, 2)
mean(silhouette(cdm_hie_org, 0.02*CDM)[, 3])
mean(silhouette(cdm_hie, sqrt(0.02*(1-0.1)*CDM**2 + 0.02*0.1*GEO**2))[, 3])
Sim(ground_truth, cdm_hie_org)
Sim(ground_truth, cdm_hie)

### LPC
choicealpha <- choicealpha(LPC, as.dist(GEO), range.alpha = alpha, K = 2, scale = F)
sil_lpc <- c()
I_lpc <- c()
for(i in 1:length(alpha)){
  D_LPC <- sqrt(0.02*(1-alpha[i])*as.matrix(LPC)**2 + 0.02*alpha[i]*GEO**2)
  lpc_tree <- hclustgeo(as.dist(LPC), as.dist(GEO), alpha = alpha[i])
  lpc_hie <- cutree(lpc_tree, 2)
  I_lpc[i] <- sum(lpc_tree$height[1:23]) 
  sil_lpc[i] <- mean(silhouette(lpc_hie, D_LPC)[, 3])
}


Q0 <- choicealpha$Q[,1]
Q1 <- choicealpha$Q[,2]
mat <- matrix(c(Q0, Q1, Q0+Q1, I_lpc), nrow = length(Q0), byrow = F)
latex_code <- xtable(mat, digits = c(0, rep(3, ncol(mat))))
print(latex_code, type = "latex")

tree <- hclustgeo(as.dist(LPC))
lpc_hie_org <- cutree(tree, 2)
lpc_tree <- hclustgeo(LPC, as.dist(GEO), alpha = 0.9)
lpc_hie <- cutree(lpc_tree, 2)
mean(silhouette(lpc_hie_org, 0.02*LPC)[, 3])
mean(silhouette(lpc_hie, sqrt(0.02*(1-0.9)*as.matrix(LPC)**2 + 0.02*0.9*GEO**2))[, 3])
Sim(ground_truth, lpc_hie_org)
Sim(ground_truth, lpc_hie)

### K = 3

choicealpha <- choicealpha(LPC, as.dist(GEO), range.alpha = alpha, K = 3, scale = F)
sil_lpc <- c()
I_lpc <- c()
for(i in 1:length(alpha)){
  D_LPC <- sqrt(0.02*(1-alpha[i])*as.matrix(LPC)**2 + 0.02*alpha[i]*GEO**2)
  lpc_tree <- hclustgeo(as.dist(LPC), as.dist(GEO), alpha = alpha[i])
  lpc_hie <- cutree(lpc_tree, 3)
  I_lpc[i] <- sum(lpc_tree$height[1:22]) 
  sil_lpc[i] <- mean(silhouette(lpc_hie, D_LPC)[, 3])
}


Q0 <- choicealpha$Q[,1]
Q1 <- choicealpha$Q[,2]
mat <- matrix(c(Q0, Q1, Q0+Q1, I_lpc), nrow = length(Q0), byrow = F)
latex_code <- xtable(mat, digits = c(0, rep(3, ncol(mat))))
print(latex_code, type = "latex")

lpc_tree <- hclustgeo(LPC, as.dist(GEO), alpha = 0.7)
lpc_hie <- cutree(lpc_tree, 3)
mean(silhouette(lpc_hie, sqrt(0.02*(1-0.7)*as.matrix(LPC)**2 + 0.02*0.7*GEO**2))[, 3])
Sim(ground_truth_new, lpc_hie)

## Ground Truth
east_coast_states <- c("CT", "DC", "DE", "FL", "MA", "MD", "ME", "NC", "NJ", 
                       "NY", "PA", "RI", "VA", "VT", "WV")
west_coast_states <- c('CA')
ground_truth <- ifelse(series_names %in% east_coast_states, 1, 2)
ground_truth_new <- c()
for(i in 1:length(series_names)){
  if(series_names[i] %in% east_coast_states){
    ground_truth_new[i] <- 1
  }else if(series_names[i] %in% west_coast_states ){
    ground_truth_new[i] <- 2
  }else{
    ground_truth_new[i] <- 3
  }
}

names(ground_truth) <- series_names
names(ground_truth_new) <- series_names
print(ground_truth)
# plot
us_states <- ne_states(country = "United States of America", returnclass = "sf")
series_names
clustering_result <- data.frame(state = series_names, cluster = lpc_hie)
state_abbrev <- data.frame(state = state.abb, name = state.name, stringsAsFactors = FALSE)
state_abbrev <- rbind(state_abbrev, data.frame(state = "DC", name = "District of Columbia"))
clustering_named <- left_join(clustering_result, state_abbrev, by = "state")
us_states_clustered <- left_join(us_states, clustering_named, by = c("name"))
us_states_clustered <- us_states_clustered %>%
  mutate(cluster = ifelse(is.na(cluster), "Unselected", as.character(cluster)))
cluster_colors <- c("1" = "#1f78b4", "2" = "#33a02c", "3" = "#e31a1c", "4" = "#ffff00","Unselected" = "gray80")

ggplot(data = us_states_clustered) +
  geom_sf(aes(fill = factor(cluster)), color = "white") +
  scale_fill_manual(values = cluster_colors, name = "Cluster")  +
  coord_sf(xlim = c(-125, -65), ylim = c(24, 50)) +
  theme_minimal() +
  labs(title = "LPC-Coord Hierarchical K = 3")
