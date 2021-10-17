
source("../R/gen_data.R")
source("../R/display.R")
source("../R/metrics.R")
library(JGL)
library(GGMselect)
library(network)
library(intergraph) #convert between network and igraph
library(BDgraph)
library(tmvtnorm)
library(matlab)
library(latex2exp)
library(viridis)
library(abind)
library(tcltk)

path <- "../figures/"
ID <- "demo"

##############################
#     generate base graph    #
#                            #
##############################
set.seed(5)
print("Base Graph")
p=20
# simulate grap
eta=0.05
base_Gr <- gen_baseGraph(p,eta)
base_gV <- network(base_Gr$G, directed=FALSE)
#convert network to igraph
ibase_gV <- asIgraph(base_gV)

plot_igraph(ibase_gV, paste(ID,"_base_graph",sep=""),path=path)
max(abs(solve(base_Gr$C)))

##############################
#       generate groups      #
#                            #
##############################
set.seed(10)
p=20
n=5000
eta=0.05
g_num = 4
par(mfrow=c(2,2), mar = c(2., 0., 0., 0.)) 
data <- list()
graph <- list()
for (i in 1:g_num){
    print(paste("Graph",i))
    
    #simulate graph and data
    group_i <- gen_groupGraph(n,p,eta,base_Gr)
    #constuct network
    gVi <- network(group_i$Gr$G, directed=FALSE)
    #convert network to igraph
    i_gVi <- asIgraph(gVi)
    
    plot_igraph(i_gVi,title=paste(ID,"_group_",i,sep=""),path=path)
    graph[i] <- list(list("igraph"=i_gVi, "ggm"=group_i$Gr))
    data[i] <-list(group_i$X)
    print(max(abs(solve(group_i$Gr$C))))
}
#png("mygraph.png")

tiff(paste(path,"demo_",ID,".tiff",sep=""), units="in", width=5, height=5, res=600)
par(mfrow=c(4,4), mar = c(2., 0., 0., 0.))
for(i in 1:4){
    plot_igraph(graph[[i]]$igraph, paste("ground truth (group ", i,")" ,sep=""),line_placement=-9)
}


print("##################################################")
print("#             Fused Graphical Lasso              #")
print("##################################################")
fgl.results=JGL(Y=data,penalty="fused",lambda1=.08,lambda2=0.05, return.whole.theta=TRUE)
print(paste("Average entropy loss:", entropy_loss(fgl.results, graph, p, g_num)))
print(paste("Average frobenious loss:", frobenious_loss(fgl.results, graph, g_num)))
print(paste("Average false positive rate:", avg_fprfnr(fgl.results, graph, g_num)[1]))
print(paste("Average false negative rate:", avg_fprfnr(fgl.results, graph, g_num)[2]))
show_result(fgl.results, graph, 4, method="fused lasso")



print("##################################################")
print("#            Group Graphical Lasso               #")
print("##################################################")
ggl.results=JGL(Y=data,penalty="group",lambda1=.08,lambda2=.09, return.whole.theta=TRUE)
print(paste("Average entropy loss:", entropy_loss(ggl.results, graph, p, g_num)))
print(paste("Average frobenious loss:", frobenious_loss(ggl.results, graph, g_num)))
print(paste("Average false positive rate:", avg_fprfnr(ggl.results, graph, g_num)[1]))
print(paste("Average false negative rate:", avg_fprfnr(ggl.results, graph, g_num)[2]))
show_result(ggl.results, graph, 4, method="group lasso")

print("##################################################")
print("#                Graphical Lasso                 #")
print("##################################################")
gl.results=JGL(Y=data,penalty="group",lambda1=.08,lambda2=.0, return.whole.theta=TRUE)
print(paste("Average entropy loss:", entropy_loss(gl.results, graph, p, g_num)))
print(paste("Average frobenious loss:", frobenious_loss(gl.results, graph, g_num)))
print(paste("Average false positive rate:", avg_fprfnr(gl.results, graph, g_num)[1]))
print(paste("Average false negative rate:", avg_fprfnr(gl.results, graph, g_num)[2]))
show_result(gl.results, graph, 4, method="lasso")

dev.off()

tiff(paste(path,"tunedAIC_",ID,".tiff",sep=""), width=5, height=5, units="in", res = 600)
par(mfrow=c(4,4), mar = c(2., 0., 0., 0.))
for(i in 1:4){
    plot_igraph(graph[[i]]$igraph, paste("ground truth group", i),line_placement=-8)
}
###########################################
#       Parameter Selection               #
#       Akaike Information Criterion      #
#       #perform 2D grid search           #
###########################################
interval_l = 10
lambda.eff <- seq(0.01, 0.3, len = interval_l)
aic_vec <- matrix(NA, length(lambda.eff),length(lambda.eff))

for(i in 1:length(lambda.eff)){
    for(j in 1:length(lambda.eff)){
        fit00 <- JGL(Y=data,penalty="fused",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
        aic_vec[i,j] <- vBIC(data, fit00, thr=0.0001)        
   }
}


tiff(file=paste(path,ID,"fgl_AIC",".tiff",sep=""), units = "in", width=5, height=5, res=600, pointsize = 12)
image(x=lambda.eff,y=lambda.eff,z=aic_vec, col=viridis(256),
      ylab=TeX('$\\lambda_2$'),xlab=TeX('$\\lambda_1$'))
contour(x=lambda.eff,y=lambda.eff,z=aic_vec, add = TRUE,nlevels = 10)
dev.off()
i_idx <- which(aic_vec == min(aic_vec), arr.ind = TRUE)[1,1]
j_idx <- which(aic_vec == min(aic_vec), arr.ind = TRUE)[1,2]

min(aic_vec)
 which(aic_vec == min(aic_vec), arr.ind = TRUE)

lam_1 <- lambda.eff[i_idx]
lam_2 <- lambda.eff[j_idx]


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
show_result(fgl_t.results, graph, 4, method="fused")

###########################################
#       Parameter Selection               #
#       Akaike Information Criterion      #
#       #perform 2D grid search           #
###########################################

aic_g_vec <- matrix(NA, length(lambda.eff),length(lambda.eff))
for(i in 1:length(lambda.eff)){
    for(j in 1:length(lambda.eff)){
        fit00 <- JGL(Y=data,penalty="group",lambda1=lambda.eff[i],lambda2=lambda.eff[j], return.whole.theta=TRUE) 
        aic_g_vec[i,j] <- vBIC(data, fit00, thr=0.0001)        
    }

}

###########################################
#         group graphical lasso           #
###########################################

tiff(file=paste(path,"ggl_AIC",".tiff",sep=""), res=600, width = 5, height = 5, units = "in", pointsize = 12)
image(x=lambda.eff, y=lambda.eff,z=aic_g_vec, col=viridis(256),
      ylab=TeX('$\\lambda_2$'),xlab=TeX('$\\lambda_1$'))
contour(x=lambda.eff,y=lambda.eff,z=aic_g_vec, add = TRUE,nlevels = 20)
dev.off()
axis(1, at = seq(0.01, 0.5, by = 0.01))
axis(2, at = seq(0.01, 0.5, by = 0.01))

#heatmap(aic_vec, scale = "row",xlab=TeX('$\\lambda_1$'))
i_idx <- which(aic_g_vec == min(aic_g_vec), arr.ind = TRUE)[1]
j_idx <- which(aic_g_vec == min(aic_g_vec), arr.ind = TRUE)[2]

min(aic_g_vec)


lam_1 <- lambda.eff[i_idx]
lam_2 <- lambda.eff[j_idx]

print(paste("lambda_1", lam_1))
print(paste("lambda_2", lam_2))



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
show_result(ggl_t.results, graph, 4, method="group")


###########################################
#       Parameter Selection               #
#       Akaike Information Criterion      #
#       #perform 2D grid search           #
###########################################

aic_gl_vec <- rep(NA, length(lambda.eff))
for(i in 1:length(lambda.eff)){
        fit00 <- JGL(Y=data,penalty="fused",lambda1=lambda.eff[i],lambda2=0., return.whole.theta=TRUE) 
        aic_gl_vec[i] <- vBIC(data, fit00, thr=0.0001)        

}

###########################################
#             graphical lasso             #
###########################################


i_idx <- which.min(aic_gl_vec) 

lam_1 <- lambda.eff[i_idx]

print(paste("optimal lambda_1", lam_1))

print("##################################################")
print("#       Graphical Lasso (tuned parameters)      #")
print("##################################################")
print(paste("optimal lambda_1", lam_1))
gl_t.results=JGL(Y=data,penalty="fused",lambda1=lam_1,lambda2=0., return.whole.theta=TRUE) 
print(paste("Average entropy loss:", entropy_loss(ggl_t.results, graph, p, 4)))
print(paste("Average frobenious loss:", frobenious_loss(gl_t.results, graph, 4)))
print(paste("Average false positive rate:", avg_fprfnr(gl_t.results, graph, 4)[1]))
print(paste("Average false negative rate:", avg_fprfnr(gl_t.results, graph, 4)[2]))
show_result(gl_t.results, graph, 4, method="gl")

dev.off()
