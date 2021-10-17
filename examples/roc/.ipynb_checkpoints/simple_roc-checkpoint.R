source("../../R/gen_data.R")
source("../../R/metrics.R")
source("../../R/display.R")
library(JGL)
library(GGMselect)
library(network)
library(intergraph) #convert between network and igraph
library(tmvtnorm)
library(latex2exp)
library(progress)

#test lambda set
lam_1 <- seq(0.01, 0.3, len = 30)
lam_2 <- c(.025, .05, .075, 1., 1.25, 1.5, 1.75, 2.0)

p <- 20
n <- 40
K <- 4
replic_num <- 100

eta_base  <- 0.2
eta_group <- 0.4


set.seed(5)
print("Base Graph")
#generate base graph
base_Gr <- gen_baseGraph(p,eta_base)

pb <- progress_bar$new(total = replic_num)
print(paste("generate", replic_num,"replications"))
graph_list <- list()
data_list  <- list()

for(q in 1:replic_num){
    set.seed(q)
    pb$tick()
    Sys.sleep(1 / replic_num)


    
    graph <- list()
    data <- list()
    for(k in 1:K){
        #simulate graph and data
        group_i <- gen_groupGraph(n, p, eta_group, base_Gr)
        #constuct network
        gVi <- network(group_i$Gr$G, directed=FALSE)
        #convert network to igraph
        i_gVi <- asIgraph(gVi)
        graph[k] <- list(list("igraph"=i_gVi, "ggm"=group_i$Gr))
        data[k] <- list(group_i$X)
    }
    data_list[q] <- list(data)
    graph_list[q] <- list(graph)
}

save(graph_list, file=paste("../../data/",replic_num,"replications_simple_graph.RData",sep=""))
save(data_list,  file=paste("../../data/",replic_num,"replications_simple_data.RData",sep=""))


load(paste("../../data/",replic_num,"replications_simple_graph.RData",sep=""))
load(paste("../../data/",replic_num,"replications_simple_data.RData",sep=""))



if(TRUE){
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
        save(fgl_roc, file=paste("../../results/FGL",replic_num,"_simple.RData",sep=""))
    }
}
}

if(TRUE){
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
            fit00 <- JGL(Y=data_list[[q]],penalty="fused",lambda1=lam_1[i],lambda2=lam_2[j], return.whole.theta=TRUE) 
            ggl_fpr[i,j] <-avg_fprfnr(fit00, graph_list[[q]], K)[1]
            ggl_fnr[i,j] <- avg_fprfnr(fit00, graph_list[[q]], K)[2]
        }
    }
    ggl_fpr_list[q] <- list(ggl_fpr)
    ggl_tpr_list[q] <- list(1-ggl_fnr)
    if(q %% 10 == 0){
        ggl_roc <- list(list("FPR"=ggl_fpr_list,"TPR"=ggl_tpr_list))
        save(ggl_roc, file=paste("../../results/GGL",replic_num,"_simple.RData",sep=""))
    }
}
}

if(TRUE){
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
        fit00 <- JGL(Y=data_list[[q]],penalty="fused",lambda1=lam_1[i],lambda2=lam_2[j], return.whole.theta=TRUE) 
        gl_fpr[i] <-avg_fprfnr(fit00, graph_list[[q]], K)[1]
        gl_fnr[i] <- avg_fprfnr(fit00, graph_list[[q]], K)[2]
    }
    gl_fpr_list[q] <- list(gl_fpr)
    gl_tpr_list[q] <- list(1-gl_fnr)
    if(q %% 10 == 0){
        gl_roc <- list(list("FPR"=gl_fpr_list,"TPR"=gl_tpr_list))
        save(gl_roc, file=paste("../../results/GL",replic_num,"_simple.RData",sep=""))
    }
}
}
