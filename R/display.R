#data visualization

plot_network <- function(network_obj){
    plot(network_obj, jitter=TRUE, usearrows = FALSE, label=1:p,displaylabels=TRUE, mode = "circle")
}


plot_igraph <- function(igraph_obj,title, path="../figures/",line_placement=-26,fname="",circ=TRUE){
    idx<-c(get.edgelist(igraph_obj))
    V(igraph_obj)$color <- "pink"
    V(igraph_obj)[idx]$color <- "skyblue"
    
    if(circ == FALSE){
        l <-  layout_with_fr
    }
    else{
        l <- layout_in_circle
    }
    
    if(fname == ""){
        fname <- title
    }

    #png(paste(path,fname,".png",sep=""),width = 480, height = 480, units = "px", pointsize = 12)
    plot(igraph_obj,vertex.color=V(igraph_obj)$color, vertex.size=15, 
         vertex.frame.color="gray", vertex.label.color="black",vertex.label.cex=0.8,layout=l)
    title(main=title,cex.main = 1, font.main= 1,line = line_placement)
    #dev.off()
}

construct_igraph <-function(jgl_matrix,n_size){
        edgeset = net.edges(jgl_matrix)[[1]]

        idx <- gsub('V','',as_ids(edgeset))
        idx_li <- unlist(strsplit(idx, "[|]"))

        gr <- make_empty_graph(n=n_size)
        gr <- gr + edges(idx_li,directed=FALSE)
        gr <-as.undirected(gr)  
        return(gr)
}

show_result <- function(est_graph, grd_graph, num, method="fused"){
#visualize result graph
    
    for (i in 1:num){
        print(paste("visualize estimated graph", i))
        
        p <- dim(est_graph$theta[[i]])[1]
        new_mat <- matrix(unlist(est_graph$theta[i]),p,p)
        adj_mat <- matrix(0,p,p)
        adj_mat[which( new_mat!= 0)] <- 1
        gr <- graph_from_adjacency_matrix(adj_mat,mode="undirected", weighted=NULL, diag=FALSE)

        #compute FN
        graph_nonzero <- length(E(grd_graph[[i]]$igraph))
        fn <-length(E(difference(grd_graph[[i]]$igraph,gr)))

        print(paste("False Negative Rate:", fn/graph_nonzero))

        #compute FP
        graph_zero <- (p+1)*p/2 - graph_nonzero
        fp <- length(E(difference(gr, grd_graph[[i]]$igraph)))
        print(paste("False Positive Rate:", fp/graph_zero))
        
        #compute Sensitivity
        #tp <- graph_nonzero -fn
        tp <- 1 - (fn/graph_nonzero)
        print(paste("Sensitivity", tp/graph_nonzero))
        
        #par(mfrow=c(1,2),mar = c(1., 0., 1., 0.),oma = c(0, 0, 0, 0))
        #plot_igraph(grd_graph[[i]]$igraph, paste("ground truth group", i))
        plot_igraph(gr, paste(method," (group ", i,")",sep=""),line_placement=-9)
        
        
    }

}

