

source("../../R/gen_data.R")
source("../../R/display.R")
source("../../R/metrics.R")
source("../../R/JGL.R")
source("../../R/admm.iters.R")
source("../../R/gete.R")
source("../../R/SSJGL.R")
source("../../R/eval.R")

library(corrplot)
library(Matrix)
 #convert between network and igraph
library(BDgraph)
library(tmvtnorm)
library(binr)

load("../../results/new_roc/roc_1/FGL20_dgp2.RData")
print("Loading FGL20_dgp2.RData...")
print(paste("Replications:", length(fgl_roc[[1]]$FPR),"lam1:",dim(fgl_roc[[1]]$FPR[[1]])[1],
            "lam2:",dim(fgl_roc[[1]]$FPR[[1]])[2]))

load("../../results/new_roc/roc_1/GGL20_dgp2.RData")
print("Loading GGL20_dgp2.RData...")
print(paste("Replications:", length(ggl_roc[[1]]$FPR),"lam1:",dim(ggl_roc[[1]]$FPR[[1]])[1],
            "lam2:",dim(ggl_roc[[1]]$FPR[[1]])[2]))

load("../../results/new_roc/roc_1/GL20_dgp2.RData")
print("Loading GL20_dgp2.RData")
print(paste("Replications:", length(gl_roc[[1]]$FPR),"lam1:",dim(gl_roc[[1]]$FPR[[1]])[1],
            "lam2:",dim(gl_roc[[1]]$FPR[[1]])[2]))

#load("../../results/scalefree_long_fused/ssjggl20_scalefree_long_fused.RData")
#print("Loading ssjggl20_scalefree_long_fused.RData...")
#print(paste("Replications:", length(ssjggl_roc[[1]]$FPR),"lam1:",dim(ssjggl_roc[[1]]$FPR[[1]])[1],
#            "lam2:",dim(ssjggl_roc[[1]]$FPR[[1]])[2]))



replic_size <- length(fgl_roc[[1]]$FPR)
lam1_size   <- dim(fgl_roc[[1]]$FPR[[1]])[1]
lam2_size   <- dim(fgl_roc[[1]]$FPR[[1]])[2]


fgl_fpr <- fgl_roc[[1]]$FPR[[1]]
fgl_tpr <- fgl_roc[[1]]$TPR[[1]]

#for development
plot(0,0,xlim = c(min(fgl_fpr),max(fgl_fpr)),ylim = c(0.,1.),type = "n")
cl <- rainbow(lam2_size)
for(l2 in 1:lam2_size){
    fpr_sort <- sort(fgl_fpr[,l2])
    tpr_sort <- fgl_tpr[,l2][order(fgl_fpr[,l2])]
    lines(fpr_sort, tpr_sort,type="l",col=cl[l2])
    
}

legend("bottomright", legend=c(0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2),
       col=c(cl), lty=1, cex=0.8)

cl <- rainbow(replic_size)
plot(0,0,xlim = c(0.,1.),ylim = c(0.,1.),type = "n")
l2 <- 1
m_fpr <- matrix(NA,replic_size, lam1_size)
m_tpr <- matrix(NA,replic_size, lam1_size)
for(s in 1:replic_size){
    fgl_fpr <- fgl_roc[[1]]$FPR[[s]]
    fgl_tpr <- fgl_roc[[1]]$TPR[[s]]
    fpr_sort <- sort(fgl_fpr[,l2])
    tpr_sort <- fgl_tpr[,l2][order(fgl_fpr[,l2])]
    lines(fpr_sort, tpr_sort,type="l",col=cl[s])
    m_fpr[s,] <- fpr_sort
    m_tpr[s,] <- tpr_sort
}

plot(0,0,xlim = c(0.,1.),ylim = c(0.,1.),type = "n")
cl <- rainbow(3)
#take average with bins 
s_fpr <- sort(m_fpr)
s_tpr <- m_tpr[order(m_fpr)]

bins_num <- 10
breaks <- seq(min(s_fpr),max(s_fpr), len=bins_num)
bins <- cut(s_fpr , breaks, include.lowest = T, right=FALSE)

b_fpr <- NULL
b_tpr <- NULL

for(i in 1:(bins_num-1)){
    if(summary(bins)[[i]]!=0){
        b_fpr <- c(b_fpr, breaks[i])
        b_tpr <- c(b_tpr, mean(s_tpr[which(findInterval(s_fpr, breaks)==i)],trim = 0.2))
    }
    
}
lines(b_fpr, b_tpr,type="l",col=cl[1])
#take average directly
d_fpr <- colMeans(m_fpr, na.rm = FALSE, dims = 1)
d_tpr2 <- apply(m_tpr, 2, mean, trim=0.2)
d_tpr1 <- colMeans(m_tpr, na.rm = FALSE, dims = 1)

lines(d_fpr, d_tpr1,type="l",col=cl[2])
lines(d_fpr, d_tpr2,type="l",col=cl[3])
legend("bottomright", legend=c("take avarage over bins","take average directly 1","take average directly 2"),
       col=c(cl), lty=1, cex=0.8)


png("../../figures/FGL_dgp2.png",width = 480, height = 480, units = "px", pointsize = 12)
plot(0,0,xlim = c(0.,1.),ylim = c(0.,1.),type = "n")

cl <- rainbow(lam2_size)

for(l2 in 1:lam2_size){
    m_fpr <- matrix(NA,replic_size, lam1_size)
    m_tpr <- matrix(NA,replic_size, lam1_size)

    for(s in 1:replic_size){

        
        fgl_fpr <- fgl_roc[[1]]$FPR[[s]]
        fgl_tpr <- fgl_roc[[1]]$TPR[[s]]
        
        fpr_sort <- sort(fgl_fpr[,l2])
        tpr_sort <- fgl_tpr[,l2][order(fgl_fpr[,l2])]
        
        #lines(fpr_sort, tpr_sort,type="l",col=cl[s])
        m_fpr[s,] <- fpr_sort
        m_tpr[s,] <- tpr_sort
        
    }
    d_fpr <- colMeans(m_fpr, na.rm = FALSE, dims = 1)
    d_tpr <- apply(m_tpr, 2, mean, trim=0.2)
    #colMeans(m_tpr, na.rm = FALSE, dims = 1)
    lines(d_fpr, d_tpr,type="l",col=cl[l2])
    #select the lam2 with largest area unterneath the ROC curve
    if(l2 == 1){
        best_fgl_fpr <- d_fpr
        best_fgl_tpr <- d_tpr 
    }

}

legend("bottomright", legend=c(0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2),
       col=c(cl), lty=1, cex=0.8)
#lines(best_fgl_fpr, best_fgl_tpr,type="l")
title(main="FGL ROC Curve",
   xlab="FPR", ylab="TPR") 
dev.off()

png("../../figures/GGL_dgp2.png",width = 480, height = 480, units = "px", pointsize = 12)
plot(NULL,xlab="FPR", ylab="TPR",xlim = c(0.,1.),ylim = c(0.,1.),type = "n")

cl <- rainbow(lam2_size)

for(l2 in 1:lam2_size){
    m_fpr <- matrix(NA,replic_size, lam1_size)
    m_tpr <- matrix(NA,replic_size, lam1_size)

    for(s in 1:replic_size){

        
        ggl_fpr <- ggl_roc[[1]]$FPR[[s]]
        ggl_tpr <- ggl_roc[[1]]$TPR[[s]]
        
        fpr_sort <- sort(ggl_fpr[,l2])
        tpr_sort <- ggl_tpr[,l2][order(ggl_fpr[,l2])]
        #lines(fpr_sort, tpr_sort,type="l",col=cl[s])
        m_fpr[s,] <- fpr_sort
        m_tpr[s,] <- tpr_sort
        
    }
    
    d_fpr <- colMeans(m_fpr, na.rm = FALSE, dims = 1)
    d_tpr <- apply(m_tpr, 2, mean, trim=0.2)
    #colMeans(m_tpr, na.rm = FALSE, dims = 1)
    lines(d_fpr, d_tpr,type="l",col=cl[l2])
    #select the lam2 with largest area unterneath the ROC curve
    if(l2 == 8){
        best_ggl_fpr <- d_fpr
        best_ggl_tpr <- d_tpr 
    }
}

legend("bottomright", legend=c(0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2),
       col=c(cl), lty=1, cex=0.8)
#lines(best_ggl_fpr, best_ggl_tpr,type="l")
title(main="GGL ROC Curve",
   xlab="FPR", ylab="TPR") 
dev.off()


plot(NULL,xlab="FPR", ylab="TPR",xlim = c(0.,1.),ylim = c(0.,1.),type = "n")

cl <- rainbow(lam2_size)

m_fpr <- matrix(NA,replic_size, 50)
m_tpr <- matrix(NA,replic_size, 50)

for(s in 1:replic_size){


    gl_fpr <- gl_roc[[1]]$FPR[[s]]
    gl_tpr <- gl_roc[[1]]$TPR[[s]]

    fpr_sort <- sort(gl_fpr)
    tpr_sort <- gl_tpr[order(gl_fpr)]

    #lines(fpr_sort, tpr_sort,type="l",col=cl[s])
    m_fpr[s,] <- fpr_sort
    m_tpr[s,] <- tpr_sort

    d_fpr <- colMeans(m_fpr, na.rm = FALSE, dims = 1)
    d_tpr <- apply(m_tpr, 2, mean, trim=0.2)
    #colMeans(m_tpr, na.rm = FALSE, dims = 1)
    lines(d_fpr, d_tpr,type="l",col=cl[1])
}

best_gl_fpr <- d_fpr
best_gl_tpr <- d_tpr



#plot(0,0,xlim = c(0.,1.),ylim = c(0.,1.),type = "n")

#cl <- rainbow(lam2_size)

#m_fpr <- matrix(NA,replic_size, lam1_size)
#m_tpr <- matrix(NA,replic_size, lam1_size)

#for(s in 1:replic_size){


#    ssjggl_fpr <- ssjggl_roc[[1]]$FPR[[s]]
#    ssjggl_tpr <- ssjggl_roc[[1]]$TPR[[s]]

#    fpr_sort <- sort(ssjggl_fpr)
#    tpr_sort <- ssjggl_tpr[order(ssjggl_fpr)]

    #lines(fpr_sort, tpr_sort,type="l",col=cl[s])
#    m_fpr[s,] <- fpr_sort
#    m_tpr[s,] <- tpr_sort

#    d_fpr <- colMeans(m_fpr, na.rm = FALSE, dims = 1)
#    d_tpr <- apply(m_tpr, 2, mean, trim=0.4)
    #colMeans(m_tpr, na.rm = FALSE, dims = 1)
#    lines(d_fpr, d_tpr,type="l",col=cl[1])
#}

#best_ssjggl_fpr <- d_fpr
#best_ssjggl_tpr <- d_tpr


png(file=paste("../../figures/dgp2_roc.png",sep=""),width = 480, height = 480, units = "px", pointsize = 12)

plot(NULL,xlab="FPR", ylab="TPR",xlim = c(0.,min(max(best_fgl_fpr),max(best_ggl_fpr),max(best_gl_fpr))),ylim = c(0.,1.),type = "n")

cl <- rainbow(3)

lines(best_fgl_fpr, best_fgl_tpr,type="l", col=cl[1])
lines(best_ggl_fpr, best_ggl_tpr,type="l", col=cl[2])
lines(best_gl_fpr, best_gl_tpr,type="l", col=cl[3])
#lines(best_ssjggl_fpr, best_ssjggl_tpr,type="l", col=cl[4])
title(main="ROC Curve") 
legend("bottomright", legend = c("FGL","GGL","GL"), col = cl, cex = 1, lwd = 2,lty=1)
dev.off()


