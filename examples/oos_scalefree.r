
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

for(ns in 1:20){
p = 100
train_num = 50
test_num = 20
K = 5
scalefree_network <- gen_scalefree_graph(train_num,p,K,rho=0.1,test=TRUE,test_num=test_num)
##############################
#       generate groups      #
#                            #
##############################

train_data <- scalefree_network[["X"]]
test_data <- scalefree_network[["testX"]]
graph <- scalefree_network[["Gr"]]


############################################
#   find optimal hyperparameters           #
#                                          #
############################################
interval_l = 20
lambda.eff <- seq(0.01, .3, len = interval_l)
fgl_bic <- matrix(NA, length(lambda.eff),length(lambda.eff))
for(i in 1:length(lambda.eff)){
    for(j in 1:length(lambda.eff)){
        fit00 <- JGL(Y=train_data,penalty="fused",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
        fgl_bic[i,j] <- vBIC(train_data, fit00)        
    }

}

i_idx <- which(fgl_bic == min(fgl_bic), arr.ind = TRUE)[1,1]
j_idx <- which(fgl_bic == min(fgl_bic), arr.ind = TRUE)[1,2]

min(fgl_bic)
 which(fgl_bic == min(fgl_bic), arr.ind = TRUE)

lam_1 <- lambda.eff[i_idx]
lam_2 <- lambda.eff[j_idx]

print(paste("lambda_1", lam_1))
print(paste("lambda_2", lam_2))


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

print(paste("lambda_1", lam_1))
print(paste("lambda_2", lam_2))

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

print(paste("lambda_1", lam_1))

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

}


