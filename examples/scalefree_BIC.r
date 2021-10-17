remove(list = ls())
source("../R/gen_data.R")
source("../R/display.R")
source("../R/metrics.R")
library(latex2exp)
library(viridis)


library(JGL)
library(BDgraph)
library(matrixcalc)
library(tmvtnorm)
library(progress)

path <- "../figures/"
ID <- "002"

n_k <- 50
K <- 5
rho <- 0.1
p <- 100

E_0 <- bdgraph.sim( p = p, graph = "scale-free" )
    
E <- graph_from_adjacency_matrix(E_0$G, mode = "undirected", weighted = NULL,
    diag = TRUE, add.colnames = NULL, add.rownames = NA)
    
#construct E_k
edge_num <- as.integer(floor(length(E(E))*rho))

graph <- list()
for(k in 1:K){
    new_E <- E_0$G
    idx <- which(new_E==0,arr.ind=TRUE)
    sample_idx <- which(idx[,1] < idx[,2])
    sample_len <- length(sample_idx)
    select <- sample(c(1:sample_len), size=edge_num, replace=FALSE)
    select_idx <- sample_idx[select]
    new_E[idx[select_idx,]] <- 1
    new_E[idx[select_idx,][,c(2,1)]] <- 1

    graph[[k]] <- list("igraph"=graph_from_adjacency_matrix(new_E, mode = "undirected", weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA))   
}
#Construct omega_k
for (k in 1:K){
    num <- length(E(graph[[k]]$igraph))
    value <- runif(num,0,1)
    value[which(value < 0.5)]  <- value[which(value < 0.5)] - 1
    new_theta <- matrix(0,p,p)


    idx <- which(as.matrix(get.adjacency(graph[[k]]$igraph))!=0,arr.ind=TRUE)
    sample_idx <- which(idx[,1] < idx[,2])
    new_theta[idx[sample_idx,]] <- value
    new_theta[idx[sample_idx,][,c(2,1)]] <- value
    diag(new_theta) <- abs(min(eigen(new_theta)$values)) + 0.1

    inv_theta <- solve(new_theta)
    diag_value <- diag(inv_theta)
    row_value <- matrix(rep(diag_value,p),p,p)
    col_value <- t(row_value)
    val <- sqrt(row_value*col_value)
    new_C <- inv_theta/val
    store <- list("G"=as.matrix(get.adjacency(graph[[k]]$igraph)),"C"=new_C,"theta"=new_theta)
    graph[[k]]$ggm <- store

}
data <- list()
for(k in 1:K){
    new_data <- rtmvnorm(n = n_k, mean = rep(0, p), sigma = graph[[k]]$ggm$C)
    #print(dim(new_data))
    data[k] <- list(new_data)
}
    

###########################################
#       Parameter Selection               #
#       Akaike Information Criterion      #
#       #perform 2D grid search           #
###########################################

print("Computing FGL BIC")
if(TRUE){
interval_l = 30
lambda.eff <- seq(0.05, 0.3, len = interval_l)
pb <- progress_bar$new(total = interval_l^2)
bic_vec <- matrix(NA, length(lambda.eff),length(lambda.eff))
for(i in 1:length(lambda.eff)){
    for(j in 1:length(lambda.eff)){
        pb$tick()
        Sys.sleep(1 / (interval_l^2))
        fit00 <- JGL(Y=data,penalty="fused",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
        bic_vec[i,j] <- vBIC(data, fit00)        
    }

}


tiff(paste(path,ID,"fgl_BIC",".tiff",sep=""), units="in", width=5, height=5, res=600, pointsize = 12)
image(x=lambda.eff,y=lambda.eff,z=bic_vec, col=viridis(256),
      ylab=TeX('$\\lambda_2$'),xlab=TeX('$\\lambda_1$'))
contour(x=lambda.eff,y=lambda.eff,z=bic_vec, add = TRUE,nlevels = 20)
dev.off()
    
i_idx <- which(bic_vec == min(bic_vec), arr.ind = TRUE)[1,1]
j_idx <- which(bic_vec == min(bic_vec), arr.ind = TRUE)[1,2]

min(bic_vec)
 which(bic_vec == min(bic_vec), arr.ind = TRUE)

lam_1 <- lambda.eff[i_idx]
lam_2 <- lambda.eff[j_idx]

print("FGL")
print(paste("lambda_1", lam_1))
print(paste("lambda_2", lam_2))

}

interval_l = 30
lambda.eff <- seq(0.05, 0.3, len = interval_l)
bic_vec <- matrix(NA, length(lambda.eff),length(lambda.eff))
pb <- progress_bar$new(total = interval_l^2)
for(i in 1:length(lambda.eff)){
    for(j in 1:length(lambda.eff)){
        pb$tick()
        Sys.sleep(1 / (interval_l^2))
        fit00 <- JGL(Y=data,penalty="group",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
        bic_vec[i,j] <- vBIC(data, fit00)        
    }

}


tiff(paste(path,ID,"ggl_BIC",".tiff",sep=""), units="in", width=5, height=5, res=600, pointsize = 12)
image(x=lambda.eff,y=lambda.eff,z=bic_vec, col=viridis(256),
      ylab=TeX('$\\lambda_2$'),xlab=TeX('$\\lambda_1$'))
contour(x=lambda.eff,y=lambda.eff,z=bic_vec, add = TRUE,nlevels = 20)

dev.off()
i_idx <- which(bic_vec == min(bic_vec), arr.ind = TRUE)[1,1]
j_idx <- which(bic_vec == min(bic_vec), arr.ind = TRUE)[1,2]

min(bic_vec)
 which(bic_vec == min(bic_vec), arr.ind = TRUE)

lam_1 <- lambda.eff[i_idx]
lam_2 <- lambda.eff[j_idx]

print("GGL")
print(paste("lambda_1", lam_1))
print(paste("lambda_2", lam_2))


interval_l = 50
lambda.eff <- seq(0.01, 0.5, len = interval_l)
pb <- progress_bar$new(total = interval_l^2)
bic_vec <- matrix(NA, length(lambda.eff))
for(i in 1:length(lambda.eff)){
        pb$tick()
        Sys.sleep(1 / (interval_l))
        fit00 <- JGL(Y=data,penalty="group",lambda1=lambda.eff[i],lambda2=0., return.whole.theta=TRUE) 
        bic_vec[i] <- vBIC(data, fit00)           
}

i_idx <- which(bic_vec == min(bic_vec))


lam_1 <- lambda.eff[i_idx]


print("GL")
print(paste("lambda_1", lam_1))
