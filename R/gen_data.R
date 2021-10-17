#generate data

gen_baseGraph <- function(p,eta){
    base_Gr <- simulateGraph(p,eta, extraeta = eta/3)
    return(base_Gr)
}

gen_groupGraph <- function(n,p,eta,base_Gr){
    Gr <- simulateGraph(p,eta, extraeta = eta/3)
    # simulate data
    Gr$G <- sign(Gr$G + base_Gr$G)
    Gr$C <- solve(solve(Gr$C) + solve(base_Gr$C)/4)
    diag(Gr$C) <- 1
    X <- rmvnorm(n, mean=rep(0,p), sigma=Gr$C)    
    return(list("Gr"=Gr, "X"=X))
}

gen_groupGraph_mask <- function(n,p,eta,base_Gr){
    
    Gr <- simulateGraph(p,eta)
    new_Gr <- list() 
    new_Gr$G <- as.matrix(bdiag(base_Gr$G, Gr$G))
    p_size <- p + dim(base_Gr$G)[1]
    
    diag(Gr$G) <- 1
    diag(base_Gr$G) <- 1
    
    #generate precision matrix
    omega_Gr <- Gr$G * solve(Gr$C)
    omega_base <- base_Gr$G * solve(base_Gr$C)
    omega <- as.matrix(bdiag(omega_Gr, omega_base))
    
    new_Gr$C <- solve(omega)
    
    # simulate data
    X <- rmvnorm(n, mean=rep(0,p_size), sigma=new_Gr$C)    
    return(list("Gr"=new_Gr, "X"=X))
}

gen_scalefree_graph <- function(n,p,K,rho=0.1,test=FALSE,test_num=0){
    E_0 <- bdgraph.sim( p = p, graph = "scale-free" )
    
    E <- graph_from_adjacency_matrix(E_0$G, mode = "undirected", weighted = NULL,
    diag = TRUE, add.colnames = NULL, add.rownames = NA)
    
    #construct E_k
    edge_num <- as.integer(floor(length(E(E))*rho))

    graph <- list()
    for(k in 1:K){
        new_E <- E_0$G
        idx <- which(new_E==0,arr.ind=TRUE)
        sample_idx <- which(idx[,1] < idx[,2])
        sample_len <- length(sample_idx)
        select <- sample(c(1:sample_len), size=edge_num, replace=FALSE)
        select_idx <- sample_idx[select]
        new_E[idx[select_idx,]] <- 1
        
        if(length(idx[select_idx,]) == 2){
            new_E[idx[select_idx,][c(2,1)]] <- 1
        }
        else{
            new_E[idx[select_idx,][,c(2,1)]] <- 1
        }

        graph[[k]] <- list("igraph"=graph_from_adjacency_matrix(new_E, mode = "undirected", weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA))   
    }
    #Construct omega_k
    for (k in 1:K){
        num <- length(E(graph[[k]]$igraph))
        value <- runif(num,0,1)
        value[which(value < 0.5)]  <- value[which(value < 0.5)] - 1
        new_theta <- matrix(0,p,p)


        idx <- which(as.matrix(get.adjacency(graph[[k]]$igraph))!=0,arr.ind=TRUE)
        sample_idx <- which(idx[,1] < idx[,2])
        new_theta[idx[sample_idx,]] <- value
        new_theta[idx[sample_idx,][,c(2,1)]] <- value
        diag(new_theta) <- abs(min(eigen(new_theta)$values)) + 0.1

        inv_theta <- solve(new_theta)
        diag_value <- diag(inv_theta)
        row_value <- matrix(rep(diag_value,p),p,p)
        col_value <- t(row_value)
        val <- sqrt(row_value*col_value)
        new_C <- inv_theta/val
        store <- list("G"=as.matrix(get.adjacency(graph[[k]]$igraph)),"C"=new_C,"theta"=new_theta)
        graph[[k]]$ggm <- store

    }
    train_data <- list()
    test_data <- list()
    for(k in 1:K){
        new_data <- rtmvnorm(n = n, mean = rep(0, p), sigma = graph[[k]]$ggm$C)
        #print(dim(new_data))
        train_data[k] <- list(new_data)

        if(test){
            new_test <- rtmvnorm(n = test_num, mean = rep(0, p), sigma = graph[[k]]$ggm$C)
            test_data[k] <- list(new_test)
        }

    }
    if(test){
        return(list("Gr"=graph, "X"=train_data, "testX"=test_data))
    }
    else{
        return(list("Gr"=graph, "X"=train_data))
    }

}

gen_scalefree_fused_graph <- function(n,p,K,rho=0.1){
    E_0 <- bdgraph.sim( p = p, graph = "scale-free" )
    
    edge_num <- as.integer(floor(sum(sum(E_0$G)/2)*rho))

    graph <- list()
    for(k in 1:K){
        new_E <- E_0$G
        idx <- which(new_E==0,arr.ind=TRUE)
        sample_idx <- which(idx[,1] < idx[,2])
        sample_len <- length(sample_idx)
        select <- sample(c(1:sample_len), size=edge_num, replace=FALSE)
        select_idx <- sample_idx[select]
        new_E[idx[select_idx,]] <- 1
        new_E[idx[select_idx,][,c(2,1)]] <- 1

        graph[[k]] <- list("igraph"=graph_from_adjacency_matrix(new_E, mode = "undirected", weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA))   
    }
    #Construct omega_k
    #construct base graph
    b_num <- as.integer(sum(E_0$G)/2)
    b_value <- runif(b_num,0,1)
    b_value[which(b_value < 0.5)]  <- b_value[which(b_value < 0.5)] - 1
    b_idx <- which(E_0$G !=0, arr.ind=TRUE)
    bs_idx <- which(b_idx[,1] < b_idx[,2])
    b_theta <- matrix(0,p,p)
    b_theta[b_idx[bs_idx,]] <- b_value
    b_theta[b_idx[bs_idx,][,c(2,1)]] <- b_value

    for (k in 1:K){
        
        g_num <- length(E(graph[[k]]$igraph)) - b_num
        value <- runif(g_num,0,1)
        value[which(value < 0.5)]  <- value[which(value < 0.5)] - 1
        new_theta <- b_theta


        idx <- which((as.matrix(get.adjacency(graph[[k]]$igraph))!=0) & (new_theta == 0), arr.ind=TRUE)
        sample_idx <- which(idx[,1] < idx[,2])
        new_theta[idx[sample_idx,]] <- value
        new_theta[idx[sample_idx,][,c(2,1)]] <- value
        
        diag(new_theta) <- abs(min(eigen(new_theta)$values)) + 0.1

        inv_theta <- solve(new_theta)
        diag_value <- diag(inv_theta)
        row_value <- matrix(rep(diag_value,p),p,p)
        col_value <- t(row_value)
        val <- sqrt(row_value*col_value)
        new_C <- inv_theta/val
        store <- list("G"=as.matrix(get.adjacency(graph[[k]]$igraph)),"C"=new_C,"theta"=new_theta)
        graph[[k]]$ggm <- store

    }
    data <- list()
    for(k in 1:K){
        new_data <- rtmvnorm(n = n, mean = rep(0, p), sigma = graph[[k]]$ggm$C)
        #print(dim(new_data))
        data[k] <- list(new_data)
    }
    return(list("Gr"=graph, "X"=data))
}

gen_random_graph <- function(n,p,K,prop=0.2,rho=0.1,test=FALSE,test_num=0){
    E_0 <- bdgraph.sim( p = p, graph = "random",prob = prop )
    
    E <- graph_from_adjacency_matrix(E_0$G, mode = "undirected", weighted = NULL,
    diag = TRUE, add.colnames = NULL, add.rownames = NA)
    
    #construct E_k
    edge_num <- as.integer(floor(length(E(E))*rho))

    graph <- list()
    for(k in 1:K){
        new_E <- E_0$G
        idx <- which(new_E==0,arr.ind=TRUE)
        sample_idx <- which(idx[,1] < idx[,2])
        sample_len <- length(sample_idx)
        select <- sample(c(1:sample_len), size=edge_num, replace=FALSE)
        select_idx <- sample_idx[select]
        new_E[idx[select_idx,]] <- 1
        
        if(length(idx[select_idx,]) == 2){
            new_E[idx[select_idx,][c(2,1)]] <- 1
        }
        else{
            new_E[idx[select_idx,][,c(2,1)]] <- 1
        }

        graph[[k]] <- list("igraph"=graph_from_adjacency_matrix(new_E, mode = "undirected", weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA))   
    }
    #Construct omega_k
    for (k in 1:K){
        num <- length(E(graph[[k]]$igraph))
        value <- runif(num,0,1)
        value[which(value < 0.5)]  <- value[which(value < 0.5)] - 1
        new_theta <- matrix(0,p,p)


        idx <- which(as.matrix(get.adjacency(graph[[k]]$igraph))!=0,arr.ind=TRUE)
        sample_idx <- which(idx[,1] < idx[,2])
        new_theta[idx[sample_idx,]] <- value
        new_theta[idx[sample_idx,][,c(2,1)]] <- value
        diag(new_theta) <- abs(min(eigen(new_theta)$values)) + 0.1

        inv_theta <- solve(new_theta)
        diag_value <- diag(inv_theta)
        row_value <- matrix(rep(diag_value,p),p,p)
        col_value <- t(row_value)
        val <- sqrt(row_value*col_value)
        new_C <- inv_theta/val
        store <- list("G"=as.matrix(get.adjacency(graph[[k]]$igraph)),"C"=new_C,"theta"=new_theta)
        graph[[k]]$ggm <- store

    }
    train_data <- list()
    test_data <- list()
    for(k in 1:K){
        new_data <- rtmvnorm(n = n, mean = rep(0, p), sigma = graph[[k]]$ggm$C)
        #print(dim(new_data))
        train_data[k] <- list(new_data)

        if(test){
            new_test <- rtmvnorm(n = test_num, mean = rep(0, p), sigma = graph[[k]]$ggm$C)
            test_data[k] <- list(new_test)
        }

    }
    if(test){
        return(list("Gr"=graph, "X"=train_data, "testX"=test_data))
    }
    else{
        return(list("Gr"=graph, "X"=train_data))
    }

}
