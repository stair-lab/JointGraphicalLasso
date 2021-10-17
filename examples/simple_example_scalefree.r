
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


path <- "../figures/simple_example_scalefree/"
ID <- "001"

n_k <- 50
K <- 5
rho <- 0.1
p <- 100

scalefree_network <- gen_scalefree_graph(n_k,p,K,rho)
graph <- scalefree_network[["Gr"]]
data  <- scalefree_network[["X"]]


if(FALSE){
print("##################################################")
print("#             Fused Graphical Lasso              #")
print("##################################################")
#lam2=0.0375
fgl.results=JGL(Y=data,penalty="fused",lambda1=.1,lambda2=.0375, return.whole.theta=TRUE)
print(paste("Average entropy loss:", entropy_loss(fgl.results, graph, p, K)))
print(paste("Average frobenious loss:", frobenious_loss(fgl.results, graph, K)))
print(paste("Average false positive rate:", avg_fprfnr(fgl.results, graph, K)[1]))
print(paste("Average false negative rate:", avg_fprfnr(fgl.results, graph, K)[2]))


print("##################################################")
print("#            Group Graphical Lasso               #")
print("##################################################")
#lam2=0.2
ggl.results=JGL(Y=data,penalty="group",lambda1=.1,lambda2=.1, return.whole.theta=TRUE)
print(paste("Average entropy loss:", entropy_loss(ggl.results, graph, p, K)))
print(paste("Average frobenious loss:", frobenious_loss(ggl.results, graph, K)))
print(paste("Average false positive rate:", avg_fprfnr(ggl.results, graph, K)[1]))
print(paste("Average false negative rate:", avg_fprfnr(ggl.results, graph, K)[2]))
#show_result(ggl.results, graph, 4, method="group")

print("##################################################")
print("#                Graphical Lasso                 #")
print("##################################################")
gl.results=JGL(Y=data,penalty="group",lambda1=.1,lambda2=.0, return.whole.theta=TRUE)
print(paste("Average entropy loss:", entropy_loss(gl.results, graph, p, K)))
print(paste("Average frobenious loss:", frobenious_loss(gl.results, graph, K)))
print(paste("Average false positive rate:", avg_fprfnr(gl.results, graph, K)[1]))
print(paste("Average false negative rate:", avg_fprfnr(gl.results, graph, K)[2]))
#show_result(gl.results, graph, 4, method="gl")
}

##############################
#     Different lambda 1     #
#     Test on all methods    #
##############################

lambda.eff <- seq(0.01, 0.3, len = 1)

el_vec <- matrix(NA,3, length(lambda.eff))
fr_vec <- matrix(NA,3, length(lambda.eff))
fpr_vec <- matrix(NA,3, length(lambda.eff))
fnr_vec <- matrix(NA,3, length(lambda.eff))
#sse_vec <- matrix(NA,3, length(lambda.eff))
for(i in 1:length(lambda.eff)){
	fit <-list()
    fit[1] <- list(JGL(Y=data,penalty="fused",lambda1=lambda.eff[i],lambda2=.2, return.whole.theta=TRUE)) 
    fit[2] <- list(JGL(Y=data,penalty="group",lambda1=lambda.eff[i],lambda2=.1, return.whole.theta=TRUE)) 
    fit[3] <- list(JGL(Y=data,penalty="group",lambda1=lambda.eff[i],lambda2=.0, return.whole.theta=TRUE)) 

    for(j in 1:3){
        fpr_vec[j,i] <-avg_fprfnr(fit[[j]], graph, K)[1]
        fnr_vec[j,i] <- avg_fprfnr(fit[[j]], graph, K)[2]
        el_vec[j,i] <- entropy_loss(fit[[j]], graph, p, K)
        fr_vec[j,i] <- frobenious_loss(fit[[j]], graph, K)
        #sse_vec[j,i] <- sum_of_squared_error(fit[[j]], graph, 4)
    }

}
#sse_vec <- sse_vec / 4

##########################
#    Entropy Loss        #
#                        #
##########################

cl <- rainbow(3)
png(paste(path,ID,"_entropy_loss",".png",sep=""),width = 480, height = 480, units = "px", pointsize = 12)

plot(0,0,xlim = c(0.0,0.3),ylim = c(0,max(el_vec)),type = "n",xlab=TeX('$\\lambda_1$'),ylab="value",main="Entropy Loss")
lines(lambda.eff, el_vec[1,], col=cl[1], type="l") 
lines(lambda.eff, el_vec[2,], col=cl[2], type="l") 
lines(lambda.eff, el_vec[3,], col=cl[3], type="l") 

legend_txt <- c("fused graphical lasso","group graphical lasso","graphical lasso")
legend("topright", legend = legend_txt, col = cl, cex = 1, lwd = 3)
dev.off()

############################
#    Frobenius Loss        #
#                          #
############################

cl <- rainbow(3)
png(paste(path,ID,"frobenius_loss",".png",sep=""),width = 480, height = 480, units = "px", pointsize = 12)

plot(0,0,xlim = c(0.0,0.3),ylim = c(0,max(fr_vec)),type = "n",xlab=TeX('$\\lambda_1$'),ylab="value",main="Frobenius Loss ")
lines(lambda.eff, fr_vec[1,], col=cl[1], type="l",lty=2) 
lines(lambda.eff, fr_vec[2,], col=cl[2], type="l",lty=2) 
lines(lambda.eff, fr_vec[3,], col=cl[3], type="l",lty=2) 
legend_txt <- c("fused graphical lasso","group graphical lasso","graphical lasso")
legend("topright", legend = legend_txt, col = cl, cex = 1, lwd = 3)
dev.off()

##########################
#          FPR,FNR       #
#                        #
##########################



cl <- rainbow(3)
png(paste(path,ID,"fpr_fnr",".png",sep=""),width = 480, height = 480, units = "px", pointsize = 12)
plot(0,0,xlim = c(0.0,0.3),ylim = c(0,max(fnr_vec)),type = "n",xlab=TeX('$\\lambda_1$'),ylab="value",main="False Positive Rate/False Negative Rate")
lines(lambda.eff, fpr_vec[1,], col=cl[1], type="l") 
lines(lambda.eff, fpr_vec[2,], col=cl[2], type="l") 
lines(lambda.eff, fpr_vec[3,], col=cl[3], type="l") 
lines(lambda.eff, fnr_vec[1,], col=cl[1], type="l",lty=2) 
lines(lambda.eff, fnr_vec[2,], col=cl[2], type="l",lty=2) 
lines(lambda.eff, fnr_vec[3,], col=cl[3], type="l",lty=2) 
legend_txt <- c("FPR - fused graphical lasso","FPR - group graphical lasso","FPR - graphical lasso",
"FNR - FPR-fused graphical lasso","FNR - group graphical lasso","FNR - graphical lasso")
legend("bottomright", legend = legend_txt, col = cl, cex = 1, lwd = 2,lty=c(1,1,1,2,2,2))

dev.off()


###########################################
#       Parameter Selection               #
#       Akaike Information Criterion      #
#       #perform 2D grid search           #
###########################################

AIC=TRUE
if(AIC){
    interval_l = 10
    lambda.eff <- seq(0.01, 0.3, len = interval_l)
    aic_vec <- matrix(NA, length(lambda.eff),length(lambda.eff))

    for(i in 1:length(lambda.eff)){
        for(j in 1:length(lambda.eff)){
            fit00 <- JGL(Y=data,penalty="fused",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
            aic_vec[i,j] <- vAIC(data, fit00)        
    }
    }

    #lambda.eff <- seq(0.01, 0.3, len = interval_l)
    #image(x=lambda.eff,y=lambda.eff,z=aic_vec, col=viridis(256),
    #      ylab=TeX('$\\lambda_2$'),xlab=TeX('$\\lambda_1$'))
    #contour(x=lambda.eff,y=lambda.eff,z=aic_vec, add = TRUE,nlevels = 10)
    #dev.copy(png, file=paste(path,ID,"fgl_AIC",".png",sep=""),width = 480, height = 480, units = "px", pointsize = 12)
    #dev.off()
    i_idx <- which(aic_vec == min(aic_vec), arr.ind = TRUE)[1,1]
    j_idx <- which(aic_vec == min(aic_vec), arr.ind = TRUE)[1,2]

    min(aic_vec)
    which(aic_vec == min(aic_vec), arr.ind = TRUE)

    lam_1 <- lambda.eff[i_idx]
    lam_2 <- lambda.eff[j_idx]

    print(paste("lambda_1", lam_1))
    print(paste("lambda_2", lam_2))
}else{
    interval_l = 10
    lambda.eff <- seq(0.01, 0.3, len = interval_l)
    bic_vec <- matrix(NA, length(lambda.eff),length(lambda.eff))
    for(i in 1:length(lambda.eff)){
        for(j in 1:length(lambda.eff)){
            fit00 <- JGL(Y=data,penalty="fused",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
            bic_vec[i,j] <- vBIC(data, fit00)        
        }

    }

    #lambda.eff <- seq(0.01, 0.3, len = interval_l)
    #image(x=lambda.eff,y=lambda.eff,z=bic_vec, col=viridis(256),
    #    ylab=TeX('$\\lambda_2$'),xlab=TeX('$\\lambda_1$'))
    #contour(x=lambda.eff,y=lambda.eff,z=bic_vec, add = TRUE,nlevels = 20)
    #dev.copy(png, file=paste(path,ID,"fgl_BIC",".png",sep=""),width = 480, height = 480, units = "px", pointsize = 12)
    #dev.off()
    i_idx <- which(bic_vec == min(bic_vec), arr.ind = TRUE)[1,1]
    j_idx <- which(bic_vec == min(bic_vec), arr.ind = TRUE)[1,2]

    min(bic_vec)
    which(bic_vec == min(bic_vec), arr.ind = TRUE)

    lam_1 <- lambda.eff[i_idx]
    lam_2 <- lambda.eff[j_idx]

    print(paste("lambda_1", lam_1))
    print(paste("lambda_2", lam_2))
        
}

#Note
#it seems like the AIC does not favor lamba_2

print("##################################################")
print("#    Fused Graphical Lasso (tuned parameters)    #")
print("##################################################")
print(paste("lambda_1", lam_1))
print(paste("lambda_2", lam_2))
fgl_t.results=JGL(Y=data,penalty="fused",lambda1=lam_1,lambda2=lam_2, return.whole.theta=TRUE) 
print(paste("Average entropy loss:", entropy_loss(fgl_t.results, graph, p, 2)))
print(paste("Average frobenious loss:", frobenious_loss(fgl_t.results, graph, 2)))
print(paste("Average false positive rate:", avg_fprfnr(fgl_t.results, graph, 2)[1]))
print(paste("Average false negative rate:", avg_fprfnr(fgl_t.results, graph, 2)[2]))
print("##################################################")
print("#             Result Visualization               #")
print("##################################################")
#show_result(fgl_t.results, graph, 4, method="fgl_tuned_001_BIC")



if(AIC){
    interval_l = 30
    lambda.eff <- seq(0.01, 0.3, len = interval_l)
    aic_g_vec <- matrix(NA, length(lambda.eff),length(lambda.eff))
    for(i in 1:length(lambda.eff)){
        for(j in 1:length(lambda.eff)){
            fit00 <- JGL(Y=data,penalty="group",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
            aic_g_vec[i,j] <- vAIC(data, fit00)        
        }

    }



    #image(x=seq(0.01,0.3,0.01),y=seq(0.01,0.3,0.01),z=aic_g_vec, col=viridis(256),
    #      ylab=TeX('$\\lambda_2$'),xlab=TeX('$\\lambda_1$'))
    #contour(x=seq(0.01,0.3,0.01),y=seq(0.01,0.3,0.01),z=aic_g_vec, add = TRUE,nlevels = 20)
    #dev.copy(png, file=paste(path,"ggl_AIC",".png",sep=""),width = 480, height = 480, units = "px", pointsize = 12)
    #dev.off()
    #axis(1, at = seq(0.01, 0.5, by = 0.01))
    #axis(2, at = seq(0.01, 0.5, by = 0.01))

    #heatmap(aic_vec, scale = "row",xlab=TeX('$\\lambda_1$'))
    i_idx <- which(aic_g_vec == min(aic_g_vec), arr.ind = TRUE)[1]
    j_idx <- which(aic_g_vec == min(aic_g_vec), arr.ind = TRUE)[2]

    min(aic_g_vec)


    lam_1 <- lambda.eff[i_idx]
    lam_2 <- lambda.eff[j_idx]

    print(paste("lambda_1", lam_1))
    print(paste("lambda_2", lam_2))
} else {
    interval_l = 30
    lambda.eff <- seq(0.01, 0.3, len = interval_l)
    bic_g_vec <- matrix(NA, length(lambda.eff),length(lambda.eff))
    for(i in 1:length(lambda.eff)){
        for(j in 1:length(lambda.eff)){
            fit00 <- JGL(Y=data,penalty="group",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
            bic_g_vec[i,j] <- vBIC(data, fit00)        
        }
    }
    #image(x=seq(0.01,0.3,0.01),y=seq(0.01,0.3,0.01),z=bic_g_vec, col=viridis(256),
    #      ylab=TeX('$\\lambda_2$'),xlab=TeX('$\\lambda_1$'))
    #contour(x=seq(0.01,0.3,0.01),y=seq(0.01,0.3,0.01),z=bic_g_vec, add = TRUE,nlevels = 20)
    #dev.copy(png, file=paste(path,"ggl_BIC",".png",sep=""),width = 480, height = 480, units = "px", pointsize = 12)
    #dev.off()
    #axis(1, at = seq(0.01, 0.5, by = 0.01))
    #axis(2, at = seq(0.01, 0.5, by = 0.01))

    #heatmap(aic_vec, scale = "row",xlab=TeX('$\\lambda_1$'))
    i_idx <- which(bic_g_vec == min(bic_g_vec), arr.ind = TRUE)[1]
    j_idx <- which(bic_g_vec == min(bic_g_vec), arr.ind = TRUE)[2]

    min(aic_g_vec)


    lam_1 <- lambda.eff[i_idx]
    lam_2 <- lambda.eff[j_idx]

    print(paste("lambda_1", lam_1))
    print(paste("lambda_2", lam_2))
}


print("##################################################")
print("#    Group Graphical Lasso (tuned parameters)    #")
print("##################################################")
print(paste("lambda_1", lam_1))
print(paste("lambda_2", lam_2))
ggl_t.results=JGL(Y=data,penalty="group",lambda1=lam_1,lambda2=lam_2, return.whole.theta=TRUE) 
print(paste("Average entropy loss:", entropy_loss(ggl_t.results, graph, p, 4)))
print(paste("Average frobenious loss:", frobenious_loss(ggl_t.results, graph, 4)))
print(paste("Average false positive rate:", avg_fprfnr(ggl_t.results, graph, 4)[1]))
print(paste("Average false negative rate:", avg_fprfnr(ggl_t.results, graph, 4)[2]))
print("##################################################")
print("#             Result Visualization               #")
print("##################################################")
#show_result(ggl_t.results, graph, 4, method="ggl_t_tuned")



interval_l = 50
lambda.eff <- seq(0.01, 1., len = interval_l)
aic_gl_vec <- rep(NA, length(lambda.eff))
for(i in 1:length(lambda.eff)){
        fit00 <- JGL(Y=data,penalty="fused",lambda1=lambda.eff[i],lambda2=0., return.whole.theta=TRUE) 
        aic_gl_vec[i] <- vAIC(data, fit00)        

}

###########################################
#             graphical lasso             #
###########################################

#image(z=aic_gl_vec, col=viridis(256),xlab=TeX('$\\lambda_1$'))
#dev.copy(png, file=paste(path,"gl_AIC",".png",sep=""),width = 480, height = 480, units = "px", pointsize = 12)
#dev.off()

i_idx <- which.min(aic_gl_vec) 

lam_1 <- lambda.eff[i_idx]

print("##################################################")
print("#       Graphical Lasso (tuned parameters)      #")
print("##################################################")
print(paste("optimal lambda_1", lam_1))
gl_t.results=JGL(Y=data,penalty="fused",lambda1=lam_1,lambda2=0., return.whole.theta=TRUE) 
print(paste("Average entropy loss:", entropy_loss(ggl_t.results, graph, p, 4)))
print(paste("Average frobenious loss:", frobenious_loss(gl_t.results, graph, 4)))
print(paste("Average false positive rate:", avg_fprfnr(gl_t.results, graph, 4)[1]))
print(paste("Average false negative rate:", avg_fprfnr(gl_t.results, graph, 4)[2]))



