remove(list = ls())
source("../../R/gen_data.R")
source("../../R/metrics.R")
source("../../R/display.R")
source("../../R/JGL.R")
source("../../R/admm.iters.R")
source("../../R/gete.R")
source("../../R/SSJGL.R")
source("../../R/eval.R")
library(BDgraph)
library(tmvtnorm)
library(JGL)
library(corrplot)
library(Matrix)
library(matrixcalc)
library(progress)

#test lambda set
lam_1 <- seq(0.0, .01, len = 10)
lam_2 <- c(0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2)
ID <- "002"

p <- 100
n_k <- 50
K <- 5
replic_num <- 20
rho <- 0.1

args <-commandArgs(trailingOnly = TRUE)

if(length(args)==0){
    stop("Please select one option: DG,FGL,GGL,GL")
}

DG = FALSE
FGL = FALSE
GGL = FALSE
GL = FALSE

if (args[1] == "DG"){
    DG = TRUE
}
if (args[1] == "FGL"){
    FGL = TRUE
}
if (args[1] == "GGL"){
    GGL = TRUE
}
if (args[1] == "GL"){
    GL = TRUE
}    
if (length(args)==2){
    p <- as.numeric(args[2])
}
pb <- progress_bar$new(total = replic_num)
print(paste("generate", replic_num,"replications"))
graph_list <- list()
data_list  <- list()

if(DG){
for(q in 1:replic_num){
    set.seed(q)
    pb$tick()
    Sys.sleep(1 / replic_num)

    E_0 <- bdgraph.sim( p = p, graph = "scale-free" )
    
    E <- graph_from_adjacency_matrix(E_0$G, mode = "undirected", weighted = NULL,
    diag = TRUE, add.colnames = NULL, add.rownames = NA)
    
    #construct E_k
    edge_num <- as.integer(floor(sum(sum(E_0$G)/2)*rho))

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
    
    data_list[q] <- list(data)
    graph_list[q] <- list(graph)
}

save(graph_list, file=paste("../../data/",replic_num,"replications_graph_dgp2.RData",sep=""))
save(data_list,  file=paste("../../data/",replic_num,"replications_data_dgp2.RData",sep=""))
print("Data saved!")
}
load(paste("../../data/",replic_num,"replications_graph_dgp2.RData",sep=""))
load(paste("../../data/",replic_num,"replications_data_dgp2.RData",sep=""))
print("Data loaded!")


if(FGL){
print("Process FGL")
fgl_fpr_list <- list() 
fgl_tpr_list <- list()

for(q in 1:replic_num){
    print(paste("processing data set",q))
    ii_size <- length(lam_1) * length(lam_2)
    fgl_fpr <- matrix(NA, length(lam_1), length(lam_2))
    fgl_fnr <- matrix(NA, length(lam_1), length(lam_2))
    pb <- progress_bar$new(total = ii_size)
    for(i in 1:length(lam_1)){
        for(j in 1:length(lam_2)){
            pb$tick()
            Sys.sleep(1 / ii_size)
            fit00 <- JGL(Y=data_list[[q]],penalty="fused",lambda1=lam_1[i],lambda2=lam_2[j], return.whole.theta=TRUE) 
            fgl_fpr[i,j] <-avg_fprfnr(fit00, graph_list[[q]], K)[1]
            fgl_fnr[i,j] <- avg_fprfnr(fit00, graph_list[[q]], K)[2]
        }
    }
    fgl_fpr_list[q] <- list(fgl_fpr)
    fgl_tpr_list[q] <- list(1-fgl_fnr)
    if(q %% 10 == 0){
        fgl_roc <- list(list("FPR"=fgl_fpr_list,"TPR"=fgl_tpr_list))
        save(fgl_roc, file=paste("../../results/FGL",replic_num,"_dgp2.RData",sep=""))
    }
}
}

if(GGL){
print("Process GGL")
ggl_fpr_list <- list() 
ggl_tpr_list <- list()

for(q in 1:replic_num){
    print(paste("processing data set",q))
    ii_size <- length(lam_1) * length(lam_2)
    ggl_fpr <- matrix(NA, length(lam_1), length(lam_2))
    ggl_fnr <- matrix(NA, length(lam_1), length(lam_2))
    pb <- progress_bar$new(total = ii_size)
    for(i in 1:length(lam_1)){
        for(j in 1:length(lam_2)){
            pb$tick()
            Sys.sleep(1 / ii_size)
            fit00 <- JGL(Y=data_list[[q]],penalty="group",lambda1=lam_1[i],lambda2=lam_2[j], return.whole.theta=TRUE) 
            ggl_fpr[i,j] <- avg_fprfnr(fit00, graph_list[[q]], K)[1]
            ggl_fnr[i,j] <- avg_fprfnr(fit00, graph_list[[q]], K)[2]
        }
    }
    ggl_fpr_list[q] <- list(ggl_fpr)
    ggl_tpr_list[q] <- list(1-ggl_fnr)
    if(q %% 10 == 0){
        ggl_roc <- list(list("FPR"=ggl_fpr_list,"TPR"=ggl_tpr_list))
        save(ggl_roc, file=paste("../../results/GGL",replic_num,"_dgp2.RData",sep=""))
    }
}
}

if(GL){
print("Process GL")
gl_fpr_list <- list() 
gl_tpr_list <- list()

for(q in 1:replic_num){
    print(paste("processing data set",q))
    ii_size <- length(lam_1) 
    gl_fpr <- matrix(NA, length(lam_1))
    gl_fnr <- matrix(NA, length(lam_1))
    pb <- progress_bar$new(total = ii_size)
    for(i in 1:length(lam_1)){

        pb$tick()
        Sys.sleep(1 / ii_size)
        fit00 <- JGL(Y=data_list[[q]],penalty="group",lambda1=lam_1[i],lambda2=0., return.whole.theta=TRUE) 
        gl_fpr[i] <-avg_fprfnr(fit00, graph_list[[q]], K)[1]
        gl_fnr[i] <- avg_fprfnr(fit00, graph_list[[q]], K)[2]
    }
    gl_fpr_list[q] <- list(gl_fpr)
    gl_tpr_list[q] <- list(1-gl_fnr)
    if(q %% 10 == 0){
        gl_roc <- list(list("FPR"=gl_fpr_list,"TPR"=gl_tpr_list))
        #save(gl_roc, file=paste("../../results/GL",replic_num,"_scalefree_",ID,".RData",sep=""))
        save(gl_roc, file=paste("../../results/GL",replic_num,"_dgp2.RData",sep=""))
    }
}
}


if(FALSE){
    print("Process SSJGGL")
    ssjggl_fpr_list <- list() 
    ssjggl_tpr_list <- list()
    for(q in 1:replic_num){
        print(paste("processing data set",q))
        lambda1 <- 1 
        lambda2 <- 1 
        v1 <- 1
        lambda.eff2 <- lambda1 + 1 + c(0:49) * 10
        v0s <- (lambda1/lambda.eff2)
        ii_size <- length(v0s)
        ssjggl_fpr <- matrix(NA, ii_size)
        ssjggl_fnr <- matrix(NA, ii_size)
        pb <- progress_bar$new(total = ii_size)
        for(i in 1:length(v0s)){
        fit00 = SSJGL(Y=data_list[[q]],penalty="group",lambda0=1, lambda1=lambda1,lambda2=lambda2, v1 = v1, v0s = v0s[i], tol.em=0.0001, a=1, b=p, doubly=TRUE, c = 0.01)
        pb$tick()
        Sys.sleep(1 / ii_size)
        ssjggl_fpr[i] <-avg_fprfnr2(fit00, graph_list[[q]], K)[1]
        ssjggl_fnr[i] <- avg_fprfnr2(fit00, graph_list[[q]], K)[2]
        }
        ssjggl_fpr_list[q] <- list(ssjggl_fpr)
        ssjggl_tpr_list[q] <- list(1-ssjggl_fnr)
        
    if(q %% 10 == 0){
        ssjggl_roc <- list(list("FPR"=ssjggl_fpr_list,"TPR"=ssjggl_tpr_list))
        save(ssjggl_roc, file=paste("../../results/ssjggl",replic_num,"_dgp2.RData",sep=""))
    }
        
    }
    
}




    
    

if(FALSE){
    print("Process SSJFGL")
    ssjfgl_fpr_list <- list() 
    ssjfgl_tpr_list <- list()
    for(q in 1:replic_num){
        print(paste("processing data set",q))
        lambda1 <- 1 
        lambda2 <- 1 
        v1 <- 1
        lambda.eff2 <- lambda1 + 1 + c(0:49) * 10
        v0s <- (lambda1/lambda.eff2)
        ii_size <- length(v0s)
        ssjfgl_fpr <- matrix(NA, ii_size)
        ssjfgl_fnr <- matrix(NA, ii_size)
        pb <- progress_bar$new(total = ii_size)
        for(i in 1:length(v0s)){
        fit00 = SSJGL(Y=data_list[[q]],penalty="fused",lambda0=1, lambda1=lambda1,lambda2=lambda2, v1 = v1, v0s = v0s[i], tol.em=0.0001, a=1, b=p, doubly=TRUE, c = 0.01)
        pb$tick()
        Sys.sleep(1 / ii_size)
        ssjfgl_fpr[i] <-avg_fprfnr2(fit00, graph_list[[q]], K)[1]
        ssjfgl_fnr[i] <- avg_fprfnr2(fit00, graph_list[[q]], K)[2]
        }
        ssjfgl_fpr_list[q] <- list(ssjfgl_fpr)
        ssjfgl_tpr_list[q] <- list(1-ssjfgl_fnr)
        
    if(q %% 10 == 0){
        ssjfgl_roc <- list(list("FPR"=ssjfgl_fpr_list,"TPR"=ssjfgl_tpr_list))
      
        save(ssjfgl_roc, file=paste("../../results/ssjfgl",replic_num,"_dgp2.RData",sep=""))
    }
        
    }
    
}
    

