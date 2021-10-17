source("../R/gen_data.R")
source("../R/display.R")
source("../R/metrics.R")
source("../R/JGL.R")
source("../R/admm.iters.R")
source("../R/gete.R")
source("../R/SSJGL.R")
source("../R/eval.R")

library(corrplot)
library(Matrix)
library(JGL)
library(GGMselect)
library(network)
library(intergraph) #convert between network and igraph
library(BDgraph)
library(tmvtnorm)
library(binr)
library(progress)

##############################
#       generate graphs      #
##############################

dataset_num <- 20


set.seed(5)
print("Base Graph")
p=20
# simulate grap
eta=0.4
base_Gr <- gen_baseGraph(p,eta)
n=40

rho <- c(0.1,0.2,0.4,0.6,0.8,1.)
label<- c(01,02,04,06,08,10)
for (r in 1:length(rho)){
    eta_g=rho[r]*eta
    if(TRUE){
    graph_list <- list()
    data_list <- list()
    
    pb <- progress_bar$new(total = dataset_num)
    print("generate dense 100 data set")
    for (k in 1:dataset_num){
        set.seed(k)
        graph <- list()
        data <- list()
        pb$tick()
        Sys.sleep(1 / dataset_num)
        for (i in 1:4){
            #simulate graph and data
            group_i <- gen_groupGraph(n,p,eta_g,base_Gr)
            #constuct network
            gVi <- network(group_i$Gr$G, directed=FALSE)
            #convert network to igraph
            i_gVi <- asIgraph(gVi)
            graph[i] <- list(list("igraph"=i_gVi, "ggm"=group_i$Gr))
            data[i] <- list(group_i$X)
        }
        ###tuning parameters and compute ROC
        graph_list[k] <- list(graph)
        data_list[k] <- list(data)
    }

    save(graph_list, file=paste("../data/20dataset_graph_rho",label[r],".RData",sep=""))
    save(data_list, file=paste("../data/20dataset_data_rho",label[r],".RData",sep=""))
    }
    load(paste("../data/20dataset_graph_rho",label[r],".RData",sep=""))
    load(paste("../data/20dataset_data_rho",label[r],".RData",sep=""))

    if(TRUE){
    print("Process FGL")
    fgl_fpr_list <- list() 
    fgl_tpr_list <- list()

    for(k in 1:dataset_num){
        print(paste("processing data set",k))
        interval_l = 30
        ii_size = interval_l^2
        lambda.eff <- seq(0.0, 0.3, len = interval_l)
        fgl_fpr <- matrix(NA, length(lambda.eff),length(lambda.eff))
        fgl_fnr <- matrix(NA, length(lambda.eff),length(lambda.eff))
        pb <- progress_bar$new(total = ii_size)
        for(i in 1:length(lambda.eff)){
            for(j in 1:length(lambda.eff)){
                pb$tick()
                Sys.sleep(1 / ii_size)
                fit00 <- JGL(Y=data_list[[k]],penalty="fused",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
                fgl_fpr[i,j] <-avg_fprfnr(fit00, graph_list[[k]], 4)[1]
                fgl_fnr[i,j] <- avg_fprfnr(fit00, graph_list[[k]], 4)[2]

            }

        }
        fgl_fpr_list[k] <- list(sort(fgl_fpr))
        fgl_tpr_list[k] <- list(1-fgl_fnr[order(fgl_fpr)])
        if(k %% 10 == 0){
            fgl_roc <- list(list("FPR"=fgl_fpr_list,"TPR"=fgl_tpr_list))
            save(fgl_roc, file=paste("../results/FGL20_rho",label[r],".RData",sep=""))
        }
    }


    print("Process GGL")
    ggl_fpr_list <- list() 
    ggl_tpr_list <- list()

    for(k in 1:dataset_num){
        print(paste("processing data set",k))
        interval_l = 30
        ii_size = interval_l^2
        lambda.eff <- seq(0.0, 0.3, len = interval_l)
        ggl_fpr <- matrix(NA, length(lambda.eff),length(lambda.eff))
        ggl_fnr <- matrix(NA, length(lambda.eff),length(lambda.eff))
        pb <- progress_bar$new(total = ii_size)
        for(i in 1:length(lambda.eff)){
            for(j in 1:length(lambda.eff)){
                pb$tick()
                Sys.sleep(1 / ii_size)
                fit00 <- JGL(Y=data_list[[k]],penalty="group",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
                ggl_fpr[i,j] <-avg_fprfnr(fit00, graph_list[[k]], 4)[1]
                ggl_fnr[i,j] <- avg_fprfnr(fit00, graph_list[[k]], 4)[2]
            }

        }
        ggl_fpr_list[k] <- list(sort(ggl_fpr))
        ggl_tpr_list[k] <- list(1-ggl_fnr[order(ggl_fpr)])
        if(k %% 10 == 0){
            ggl_roc <- list(list("FPR"=ggl_fpr_list,"TPR"=ggl_tpr_list))
            save(ggl_roc, file=paste("../results/GGL20_rho",label[r],".RData",sep=""))
        }
    }
    }

    print("Process GL")
    gl_fpr_list <- list() 
    gl_tpr_list <- list()

    for(k in 1:dataset_num){
        print(paste("processing data set",k))
        interval_l = 100
        ii_size = interval_l
        lambda.eff <- seq(0.0, 0.5, len = interval_l)
        gl_fpr <- matrix(NA, length(lambda.eff))
        gl_fnr <- matrix(NA, length(lambda.eff))
        pb <- progress_bar$new(total = ii_size)
        for(i in 1:length(lambda.eff)){
                pb$tick()
                Sys.sleep(1 / ii_size)
                fit00 <- JGL(Y=data_list[[k]],penalty="group",lambda1=lambda.eff[i],lambda2=0., return.whole.theta=TRUE) 
                gl_fpr[i] <-avg_fprfnr(fit00, graph_list[[k]], 4)[1]
                gl_fnr[i] <- avg_fprfnr(fit00, graph_list[[k]], 4)[2]

        }
        gl_fpr_list[k] <- list(sort(gl_fpr))
        gl_tpr_list[k] <- list(1-gl_fnr[order(gl_fpr)])
        if(k %% 10 == 0){
            gl_roc <- list(list("FPR"=gl_fpr_list,"TPR"=gl_tpr_list))
            save(gl_roc, file=paste("../results/GL20_rho",label[r],".RData",sep=""))
        }
    }
}