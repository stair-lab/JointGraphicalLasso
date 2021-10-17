remove(list = ls())
source("../R/JGL.R")
source("../R/admm.iters.R")
source("../R/gete.R")
source("../R/SSJGL.R")
source("../R/eval.R")
source("../R/metrics.R")
library(latex2exp)
library(BDgraph)
library(tmvtnorm)
library(JGL)
library(corrplot)
library(Matrix)
library(network)
library(intergraph)
library(progress)

###################################################
##  Simulate data
###################################################



set.seed(1)

graph_list <- list()
data_list <- list()
dataset_num <- 100
pb <- progress_bar$new(total = dataset_num)
print(paste("generate",dataset_num,"data set"))


for(k in 1:dataset_num){
    set.seed(k)
    pb$tick()
    Sys.sleep(1 / dataset_num)
    P <- 10 #dimension
    PP <- 100
    N <- rep(150, 2)

    rho <- runif(2, min = 0, max = 1)
    remove <- 4
    data <- g <- thetatrue <- NULL

    #generate s1, #G is the adjacent matirx
    s1 <- bdgraph.sim.rho(n = 1, p = P, vis = FALSE, graph = "AR1")$G #vis:visualize the true graph
    diag(s1) <- 0

    #generate s2
    s2 <- bdgraph.sim.rho(n = 1, p = P/2, vis = FALSE, graph = "AR1")$G
    diag(s2) <- 0
    s2 <- bdiag(s2, diag(0, P/2))

    s0 <- matrix(0, PP - P, PP - P)
    #bdiag means block diagonal matrix
    g[[1]] <- as.matrix(bdiag(s1, s0))
    g[[2]] <- as.matrix(bdiag(s2, s0))

    diag(s1) <- diag(s2) <- 1

    #K is the precision matrix
    theta1 <- s1 * bdgraph.sim.rho(n = 1, p = P, vis = FALSE, graph = "AR1", rho=rho[1])$K
    theta2 <- s2 * bdgraph.sim.rho(n = 1, p = P, vis = FALSE, graph = "AR1", rho=rho[2])$K
    thetatrue[[1]] <- (as.matrix(bdiag(theta1, diag(1, PP-P))))
    thetatrue[[2]] <- (as.matrix(bdiag(theta2, diag(1, PP-P))))

    graph <- list()

    tmp_1 <- list("C"=solve(thetatrue[[1]]),"G"=sign(abs(thetatrue[[1]])))
    gV1 <- network(tmp_1$G, directed=FALSE)
    #convert network to igraph
    i_gV1 <- asIgraph(gV1)
    graph[[1]] <- list("ggm"=tmp_1,"igraph"=i_gV1) 

    tmp_2 <- list("C"=solve(thetatrue[[2]]),"G"=sign(abs(thetatrue[[2]])))
    gV2 <- network(tmp_2$G, directed=FALSE)
    #convert network to igraph
    i_gV2 <- asIgraph(gV2)
    graph[[2]] <- list("ggm"=tmp_2,"igraph"=i_gV2) 


    #solve(..), compute the inverse
    data[[1]] <- rtmvnorm(n = N[1], mean = rep(0, PP), sigma = solve(thetatrue[[1]]))
    data[[2]] <- rtmvnorm(n = N[2], mean  = rep(0, PP), sigma = solve(thetatrue[[2]]))
    data_list[k] <- list(data)
    graph_list[k] <- list(graph)

}

save(graph_list, file=paste("../data/",dataset_num,"dataset_graph_ssjgl.RData", sep=""))
save(data_list,  file=paste("../data/",dataset_num,"dataset_data_ssjgl.RData",sep=""))

#rm(list=setdiff(ls(), "dataset_num"))

load(paste("../data/",dataset_num,"dataset_graph_ssjgl.RData", sep=""))
load(paste("../data/",dataset_num,"dataset_data_ssjgl.RData",sep=""))


if(TRUE){
print("Process FGL")
fgl_fpr_list <- list() 
fgl_tpr_list <- list()

for(k in 1:dataset_num){
    print(paste("processing data set",k))
    interval_l = 30
    ii_size = interval_l^2
    lambda.eff <- seq(0.01, 0.5, len = interval_l)
    fgl_fpr <- matrix(NA, length(lambda.eff),length(lambda.eff))
    fgl_fnr <- matrix(NA, length(lambda.eff),length(lambda.eff))
    pb <- progress_bar$new(total = ii_size)
    for(i in 1:length(lambda.eff)){
        for(j in 1:length(lambda.eff)){
            pb$tick()
            Sys.sleep(1 / ii_size)
            fit00 <- JGL(Y=data_list[[k]],penalty="fused",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
            fgl_fpr[i,j] <-avg_fprfnr(fit00, graph_list[[k]], 2)[1]
            fgl_fnr[i,j] <- avg_fprfnr(fit00, graph_list[[k]], 2)[2]
        }
    }
    fgl_fpr_list[k] <- list(sort(fgl_fpr))
    fgl_tpr_list[k] <- list(1-fgl_fnr[order(fgl_fpr)])
    if(k %% 10 == 0){
        fgl_roc <- list(list("FPR"=fgl_fpr_list,"TPR"=fgl_tpr_list))
        save(fgl_roc, file=paste("../results/FGL",dataset_num,"_ssjgl.RData",sep=""))
    }
}
}

if(TRUE){
print("Process GGL")
ggl_fpr_list <- list() 
ggl_tpr_list <- list()

for(k in 1:dataset_num){
    print(paste("processing data set",k))
    interval_l = 30
    ii_size = interval_l^2
    lambda.eff <- seq(0.01, 0.5, len = interval_l)
    ggl_fpr <- matrix(NA, length(lambda.eff),length(lambda.eff))
    ggl_fnr <- matrix(NA, length(lambda.eff),length(lambda.eff))
    pb <- progress_bar$new(total = ii_size)
    for(i in 1:length(lambda.eff)){
        for(j in 1:length(lambda.eff)){
            pb$tick()
            Sys.sleep(1 / ii_size)
            fit00 <- JGL(Y=data_list[[k]],penalty="group",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
            ggl_fpr[i,j] <-avg_fprfnr(fit00, graph_list[[k]], 2)[1]
            ggl_fnr[i,j] <- avg_fprfnr(fit00, graph_list[[k]], 2)[2]
        }
    }
    ggl_fpr_list[k] <- list(sort(ggl_fpr))
    ggl_tpr_list[k] <- list(1-ggl_fnr[order(ggl_fpr)])
    if(k %% 10 == 0){
        ggl_roc <- list(list("FPR"=ggl_fpr_list,"TPR"=ggl_tpr_list))
        save(ggl_roc, file=paste("../results/GGL",dataset_num,"_ssjgl.RData",sep=""))
    }
}}


if(TRUE){
print("Process GL")
gl_fpr_list <- list() 
gl_tpr_list <- list()

for(k in 1:dataset_num){
    print(paste("processing data set",k))
    interval_l = 200
    ii_size = interval_l
    lambda.eff <- seq(0.0, 0.5, len = interval_l)
    gl_fpr <- matrix(NA, length(lambda.eff))
    gl_fnr <- matrix(NA, length(lambda.eff))
    pb <- progress_bar$new(total = ii_size)
    for(i in 1:length(lambda.eff)){
            pb$tick()
            Sys.sleep(1 / ii_size)
            fit00 <- JGL(Y=data_list[[k]],penalty="group",lambda1=lambda.eff[i],lambda2=0., return.whole.theta=TRUE) 
            gl_fpr[i] <-avg_fprfnr(fit00, graph_list[[k]], 2)[1]
            gl_fnr[i] <- avg_fprfnr(fit00, graph_list[[k]], 2)[2]
    }
    gl_fpr_list[k] <- list(sort(gl_fpr))
    gl_tpr_list[k] <- list(1-gl_fnr[order(gl_fpr)])
    if(k %% 10 == 0){
        gl_roc <- list(list("FPR"=gl_fpr_list,"TPR"=gl_tpr_list))
        save(gl_roc, file=paste("../results/GL",dataset_num,"_ssjgl.RData",sep=""))
    }
}}
