
source("../R/gen_data.R")
source("../R/display.R")
source("../R/metrics.R")
library(JGL)
library(GGMselect)
library(network)
library(intergraph) #convert between network and igraph
library(BDgraph)
library(tmvtnorm)
library(viridis)
library(latex2exp)
library(fields)

############################################
#   compute out of sample loglikelihood    #
#                                          #
############################################

oos_loglikelihood <- function(x_test, est_graph){
    num <- length(x_test)
    oos_sum <- 0.
    for(i in 1:num){
        sample_cov <- cov(x_test[[i]], x_test[[i]])
        log_det <- determinant(est_graph$theta[[i]], logarithm = TRUE)$modulus[1][1]
        tr_sum <- sum(diag(sample_cov %*% est_graph$theta[[i]]))
        
        oos_sum <- oos_sum + (log_det - tr_sum)
    }
    return(oos_sum)
}

gen_groupGraph_with_test <- function(train_num,test_num,p,eta,base_Gr){
    Gr <- simulateGraph(p,eta, extraeta = eta/3)
    # simulate data
    Gr$C <- solve(solve(Gr$C) + solve(base_Gr$C)/4)
    diag(Gr$C) <- 1
    X_train <- rmvnorm(train_num, mean=rep(0,p), sigma=Gr$C)
    X_test <- rmvnorm(test_num, mean=rep(0,p), sigma=Gr$C)
    Gr$G <- sign(Gr$G + base_Gr$G)
    return(list("Gr"=Gr, "X_train"=X_train, "X_test"=X_test))
}

p = 20
train_num = 40
test_num = 20

set.seed(5)
print("Base Graph")
eta=0.2
base_Gr <- gen_baseGraph(p,eta)
base_gV <- network(base_Gr$G, directed=FALSE)
#convert network to igraph
ibase_gV <- asIgraph(base_gV)
#plot_igraph(ibase_gV, "base graph", save=TRUE,path="../figures/model_selection")

##############################
#       generate groups      #
#                            #
##############################
set.seed(10)

eta=0.1
par(mfrow=c(2,2), mar = c(2., 0., 0., 0.)) 
train_data <- list()
test_data <- list()
graph <- list()
for (i in 1:4){
    print(paste("Graph",i))
    
    #simulate graph and data
    group_i <- gen_groupGraph_with_test (train_num,test_num,p,eta,base_Gr)
    
    #constuct network
    gVi <- network(group_i$Gr$G, directed=FALSE)
    
    #convert network to igraph
    i_gVi <- asIgraph(gVi)
    
    #plot_igraph(i_gVi,title=paste("group",i),save=TRUE,line_placement=-19)
    graph[i] <- list(list("igraph"=i_gVi, "ggm"=group_i$Gr))
    train_data[i] <-list(group_i$X_train)
    test_data[i] <- list(group_i$X_test)
}


############################################
#   find optimal hyperparameters           #
#                                          #
############################################
interval_l = 50
lambda.eff <- seq(0.01, 0.5, len = interval_l)
fgl_bic <- matrix(NA, length(lambda.eff),length(lambda.eff))
for(i in 1:length(lambda.eff)){
    for(j in 1:length(lambda.eff)){
        fit00 <- JGL(Y=train_data,penalty="fused",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
        fgl_bic[i,j] <- vBIC(train_data, fit00)        
    }

}


i_idx <- which(fgl_bic == min(fgl_bic), arr.ind = TRUE)[1,1]
j_idx <- which(fgl_bic == min(fgl_bic), arr.ind = TRUE)[1,2]


lam_1 <- lambda.eff[i_idx]
lam_2 <- lambda.eff[j_idx]


print("##################################################\n")
print("#             Fused Graphical Lasso              #\n")
print("##################################################\n")
print(paste("lambda_1", lam_1,"\n"))
print(paste("lambda_2", lam_2,"\n"))
fgl.results=JGL(Y=train_data,penalty="fused",lambda1=lam_1,lambda2=lam_2, return.whole.theta=TRUE)
print(paste("Average entropy loss:", entropy_loss(fgl.results, graph, p, 4),"\n"))
print(paste("Average frobenious loss:", frobenious_loss(fgl.results, graph, 4),"\n"))
print(paste("Average false positive rate:", avg_fprfnr(fgl.results, graph, 4)[1],"\n"))
print(paste("Average false negative rate:", avg_fprfnr(fgl.results, graph, 4)[2],"\n"))


print(paste('FGL out of sample likelihood:', oos_loglikelihood(test_data,fgl.results),"\n"))

############################################
#   find optimal hyperparameters           #
#                                          #
############################################


ggl_bic <- matrix(NA, length(lambda.eff),length(lambda.eff))
for(i in 1:length(lambda.eff)){
    for(j in 1:length(lambda.eff)){
        fit00 <- JGL(Y=train_data,penalty="group",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
        ggl_bic[i,j] <- vBIC(train_data, fit00)        
    }

}


#heatmap(bic_vec, scale = "row",xlab=TeX('$\\lambda_1$'))
i_idx <- which(ggl_bic == min(ggl_bic), arr.ind = TRUE)[1,1]
j_idx <- which(ggl_bic == min(ggl_bic), arr.ind = TRUE)[1,2]

#min(fgl_bic)==fgl_bic[i_idx,j_idx]

which(ggl_bic == min(ggl_bic), arr.ind = TRUE)
lam_1 <- lambda.eff[i_idx]
lam_2 <- lambda.eff[j_idx]


print("##################################################\n")
print("#             Group Graphical Lasso              #\n")
print("##################################################\n")
print(paste("lambda_1", lam_1,"\n"))
print(paste("lambda_2", lam_2,"\n"))
ggl.results=JGL(Y=train_data,penalty="group",lambda1=lam_1,lambda2=lam_2, return.whole.theta=TRUE)
print(paste("Average entropy loss:", entropy_loss(ggl.results, graph, p, 4),"\n"))
print(paste("Average frobenious loss:", frobenious_loss(ggl.results, graph, 4),"\n"))
print(paste("Average false positive rate:", avg_fprfnr(ggl.results, graph, 4)[1],"\n"))
print(paste("Average false negative rate:", avg_fprfnr(ggl.results, graph, 4)[2],"\n"))

print(paste("GGL out of sample likelihood:",oos_loglikelihood(test_data,ggl.results),"\n"))

############################################
#   find optimal hyperparameters           #
#                                          #
############################################



gl_bic <- matrix(NA, length(lambda.eff))
for(i in 1:length(lambda.eff)){
    fit00 <- JGL(Y=train_data,penalty="group",lambda1=lambda.eff[i],lambda2=0., return.whole.theta=TRUE) 
    gl_bic[i] <- vBIC(train_data, fit00)        
}


i_idx <- which(gl_bic == min(gl_bic), arr.ind = TRUE)[1]

lam_1 <- lambda.eff[i_idx]


#########################################
#          test graphical lasso         #
#                                       #
#########################################

print("##################################################\n")
print("#                Graphical Lasso                 #\n")
print("##################################################\n")
print(paste("lambda_1", lam_1,"\n"))
gl.results=JGL(Y=train_data,penalty="fused",lambda1=.08,lambda2=.0, return.whole.theta=TRUE)
print(paste("Average entropy loss:", entropy_loss(gl.results, graph, p, 4),"\n"))
print(paste("Average frobenious loss:", frobenious_loss(gl.results, graph, 4),"\n"))
print(paste("Average false positive rate:", avg_fprfnr(gl.results, graph, 4)[1],"\n"))
print(paste("Average false negative rate:", avg_fprfnr(gl.results, graph, 4)[2],"\n"))


print(paste("GL out of sample likelihood:",oos_loglikelihood(test_data,gl.results),"\n"))




