
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

#todo list test on 20 replications
path <- "../figures/"
ID <- "005Small"
n_k <- 10
K <- 5
rho <- 0.3
p <- c(20, 50, 100, 150, 200)

entropy_metric <- matrix(0, 3, length(p))
frobenious_metric <- matrix(0, 3, length(p))
fnr_metric <- matrix(0, 3, length(p))
fpr_metric <- matrix(0, 3, length(p))

if(FALSE){
for(r in 1:length(p)){
    print(p[r])
    scalefree_network <- gen_scalefree_graph(n_k,p[r],K,rho)
    graph <- scalefree_network[["Gr"]]
    data  <- scalefree_network[["X"]]

    interval_l = 20
    lambda1.eff <- seq(0.01, 0.5, len = interval_l)
    lambda2.eff <- c(0.05, 0.1, 0.15, 0.2)
    bic_fgl <- matrix(NA, length(lambda1.eff),length(lambda2.eff))
    bic_ggl <- matrix(NA, length(lambda1.eff),length(lambda2.eff))
    bic_gl <-  matrix(NA, length(lambda1.eff))
    itr <- length(lambda1.eff)*length(lambda2.eff)
    pb <- progress_bar$new(total = itr)
    for(i in 1:length(lambda1.eff)){
        for(j in 1:length(lambda2.eff)){
            pb$tick()
            Sys.sleep(1 / itr)
            fit_fused <- JGL(Y=data,penalty="fused",lambda1=lambda1.eff[i],lambda2=lambda2.eff[j], return.whole.theta=TRUE) 
            fit_group <- JGL(Y=data,penalty="group",lambda1=lambda1.eff[i],lambda2=lambda2.eff[j], return.whole.theta=TRUE) 
            bic_fgl[i,j] <- vBIC(data, fit_fused)
            bic_ggl[i,j] <- vBIC(data, fit_group)
        }
        fit_gl <- JGL(Y=data,penalty="group",lambda1=lambda1.eff[i],lambda2=0., return.whole.theta=TRUE) 
        bic_gl[i] <- vBIC(data, fit_gl)
    }
    #the coefficients for fgl
    i_idx <- which(bic_fgl == min(bic_fgl), arr.ind = TRUE)[1,1]
    j_idx <- which(bic_fgl == min(bic_fgl), arr.ind = TRUE)[1,2]

    fgl_1 <- lambda1.eff[i_idx]
    fgl_2 <- lambda2.eff[j_idx]

    #the coefficients for ggl
    i_idx <- which(bic_ggl == min(bic_ggl), arr.ind = TRUE)[1,1]
    j_idx <- which(bic_ggl == min(bic_ggl), arr.ind = TRUE)[1,2]

    ggl_1 <- lambda1.eff[i_idx]
    ggl_2 <- lambda2.eff[j_idx]

    #the coefficients fo gl
    i_idx <- which.min(bic_gl) 

    gl_1 <- lambda1.eff[i_idx]

    #test on optimal model
    fit <-list()
    fit[1] <- list(JGL(Y=data, penalty="fused", lambda1=fgl_1,lambda2=fgl_2, return.whole.theta=TRUE)) 
    fit[2] <- list(JGL(Y=data, penalty="group", lambda1=ggl_1,lambda2=ggl_2, return.whole.theta=TRUE)) 
    fit[3] <- list(JGL(Y=data, penalty="group", lambda1=gl_1, lambda2=0.,    return.whole.theta=TRUE)) 
    
    
    for(j in 1:3){
        fpr_metric[j,r] <-avg_fprfnr(fit[[j]], graph, K)[1]
        fnr_metric[j,r] <- avg_fprfnr(fit[[j]], graph, K)[2]
        entropy_metric[j,r] <- entropy_loss(fit[[j]], graph, p[r], K)
        frobenious_metric[j,r] <- frobenious_loss(fit[[j]], graph, K)
    }
}

metrics <- list("FPR"=fpr_metric,"FNR"=fnr_metric,"FL"=frobenious_metric,"EL"=entropy_metric)
save(metrics, file=paste("../results/loss_graphsize_scalefree/",ID,".RData",sep=""))

}

load(paste("../results/loss_graphsize_scalefree/",ID,".RData",sep=""))

cl <- rainbow(3)
tiff(paste(path,ID,"_entropy_loss",".tiff",sep=""), units="in", width=5, height=5, res=600, pointsize = 12)

plot(0,0,xlim=c(min(log10(p)),max(log10(p))), ylim = c(0,max(metrics[["EL"]])),type = "n",xlab=TeX('log graph size'),ylab="value",main="Entropy Loss")
lines(log10(p), metrics[["EL"]][1,], col=cl[1], type="l") 
lines(log10(p), metrics[["EL"]][2,], col=cl[2], type="l") 
lines(log10(p), metrics[["EL"]][3,], col=cl[3], type="l") 

legend_txt <- c("fused graphical lasso","group graphical lasso","graphical lasso")
legend("topright", legend = legend_txt, col = cl, cex = 1, lwd = 3)
dev.off()



png(paste(path,ID,"frobenius_loss",".png",sep=""),width = 480, height = 480, units = "px", pointsize = 12)

plot(0,0, xlim=c(min(log10(p)),max(log10(p))), ylim = c(0,max(metrics[["FL"]][1,])),type = "n",xlab=TeX('log graph size'),ylab="value",main="Frobenius Loss ")
lines(log10(p), metrics[["FL"]][1,], col=cl[1], type="l",lty=2) 
lines(log10(p), metrics[["FL"]][2,], col=cl[2], type="l",lty=2) 
lines(log10(p), metrics[["FL"]][3,], col=cl[3], type="l",lty=2) 
legend_txt <- c("fused graphical lasso","group graphical lasso","graphical lasso")
legend("topright", legend = legend_txt, col = cl, cex = 1, lwd = 3)
dev.off()


png(paste(path,ID,"fpr_fnr",".png",sep=""),width = 480, height = 480, units = "px", pointsize = 12)
plot(0,0, xlim=c(min(log10(p)),max(log10(p))), ylim = c(0.,1.),type = "n",xlab=TeX('log graph size'),ylab="value",main="False Positive Rate/False Negative Rate")
lines(log10(p), metrics[["FPR"]][1,], col=cl[1], type="l") 
lines(log10(p), metrics[["FPR"]][2,], col=cl[2], type="l") 
lines(log10(p), metrics[["FPR"]][3,], col=cl[3], type="l") 
lines(log10(p), metrics[["FNR"]][1,], col=cl[1], type="l",lty=2) 
lines(log10(p), metrics[["FNR"]][2,], col=cl[2], type="l",lty=2) 
lines(log10(p), metrics[["FNR"]][3,], col=cl[3], type="l",lty=2) 
legend_txt <- c("FPR - fused graphical lasso","FPR - group graphical lasso","FPR - graphical lasso",
"FNR - FPR-fused graphical lasso","FNR - group graphical lasso","FNR - graphical lasso")
legend("bottomright", legend = legend_txt, col = cl, cex = 1, lwd = 2,lty=c(1,1,1,2,2,2))

dev.off()
