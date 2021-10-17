#evaluation 
entropy_loss <- function(est_graph, grd_graph, p, num ){
    el_sum <- 0.
    for(i in 1:num){
        tmp_matrix <- grd_graph[[i]]$ggm$C %*% est_graph$theta[[i]]
        tmp_tr <- sum(diag(tmp_matrix))
        tmp_det <- determinant(tmp_matrix, logarithm = TRUE)
        el_sum <- el_sum + tmp_tr - tmp_det$modulus[1][1]
    }
    el_sum <- el_sum / num
    el_sum <- el_sum - p
    return(el_sum)
}

frobenious_loss <- function(est_graph, grd_graph, num){
    fl_sum <- 0.
    for (i in 1:num){
        grd_c <- solve(grd_graph[[i]]$ggm$C)
        sub_matrix <- grd_c - est_graph$theta[[i]]
        sub_norm2 <- norm(sub_matrix, "F") ^ 2
        grd_norm2 <- norm(grd_c, "F") ^ 2
        fl_sum <- fl_sum + (sub_norm2/grd_norm2)
    }
    fl_sum <- fl_sum / num
    return(fl_sum)
}
sum_of_squared_error <- function(est_graph, grd_graph, num){
    fl_sum <- 0.
    for (i in 1:num){
        grd_c <- solve(grd_graph[[i]]$ggm$C)
        sub_matrix <- grd_c - est_graph$theta[[i]]
        diag(sub_matrix)<- 0.
        sub_norm2 <- norm(sub_matrix, "F") ^ 2
        fl_sum <- fl_sum + sub_norm2
    }
    return(fl_sum)
}
avg_fprfnr <- function(est_graph, grd_graph, num, thr=0.001){
    fn_sum <- 0.
    fp_sum <- 0.
    for(i in 1:num){
        p <- dim(est_graph$theta[[i]])[1]
        new_mat <- matrix(unlist(est_graph$theta[[i]]),p,p)
        adj_mat <- matrix(0,p,p)
        adj_mat[which( abs(new_mat) >= thr)] <- 1
        gr <- graph_from_adjacency_matrix(adj_mat,mode="undirected", weighted=NULL, diag=FALSE)
        #gr <- construct_igraph(est_graph$theta[i])
        
        graph_nonzero <- length(E(grd_graph[[i]]$igraph))
        fn <-length(E(difference(grd_graph[[i]]$igraph,gr)))
        fn_sum <- fn_sum + fn/graph_nonzero

        graph_zero <- (p+1)*p/2 - graph_nonzero
        fp <- length(E(difference(gr, grd_graph[[i]]$igraph)))
        fp_sum <- fp_sum + fp/graph_zero
    }
    fn_sum <- fn_sum / num
    fp_sum <- fp_sum / num
    return(c(fp_sum,fn_sum))
}
avg_fprfnr2 <- function(est_graph, grd_graph, num){
    fn_sum <- 0.
    fp_sum <- 0.
    for(i in 1:num){
        p <- dim(est_graph$theta[[1]][[i]])[1]
        new_mat <- matrix(unlist(est_graph$theta[[1]][i]),p,p)
        adj_mat <- matrix(0,p,p)
        adj_mat[which( new_mat!= 0)] <- 1
        gr <- graph_from_adjacency_matrix(adj_mat,mode="undirected", weighted=NULL, diag=FALSE)
        #gr <- construct_igraph(est_graph$theta[i])
        
        graph_nonzero <- length(E(grd_graph[[i]]$igraph))
        fn <-length(E(difference(grd_graph[[i]]$igraph,gr)))
        fn_sum <- fn_sum + fn/graph_nonzero

        graph_zero <- (p+1)*p/2 - graph_nonzero
        fp <- length(E(difference(gr, grd_graph[[i]]$igraph)))
        fp_sum <- fp_sum + fp/graph_zero
    }
    fn_sum <- fn_sum / num
    fp_sum <- fp_sum / num
    return(c(fp_sum,fn_sum))
}

vAIC <- function(X, est_graph,thr=0.001){
    num <- length(X)
    AIC_acc <- 0.
    for(i in 1:num){
        
        data_num <- dim(X[[i]])[1]
        sample_cov <- cov(X[[i]], X[[i]])
        tr_sum <- sum(diag(sample_cov %*% est_graph$theta[[i]]))

        log_det <- determinant(est_graph$theta[[i]], logarithm = TRUE)$modulus[1][1]

        E <- sum(sum(abs(est_graph$theta[[i]]) >= thr))
        AIC_acc <- AIC_acc + data_num*(tr_sum - log_det) + 2*E
    }
    return(AIC_acc)
}
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
vBIC <- function(X, est_graph, thr=0.001){
    num <- length(X)
    BIC_acc <- 0.
    for(i in 1:num){
        
        data_num <- dim(X[[i]])[1]
        sample_cov <- cov(X[[i]], X[[i]])
        tr_sum <- sum(diag(sample_cov %*% est_graph$theta[[i]]))

        log_det <- determinant(est_graph$theta[[i]], logarithm = TRUE)$modulus[1][1]

        E <- sum(sum(abs(est_graph$theta[[i]]) >= thr))
        BIC_acc <- BIC_acc + (tr_sum - log_det) + (log(data_num)*E/data_num)
    }
    return(BIC_acc)
}

###########################################
#       Parameter Selection               #
#       Cross Validation                  #
###########################################
convert2list <- function(X,len){
    new_list <- list()
    for (i in 1:len){
        new_matrix <- X[,,i]
        new_list[i] <- list(new_matrix)
    }
    return(new_list)
}
CV <- function(X, D, method="fused", lam1, lam2){
    num <- length(X)
    nn <- dim(X[[1]])[1]
    pp <- dim(X[[1]])[2]
    CV_acc <- 0.
    X_matrix <- array(unlist(X),dim=c(nn,pp,length(X)))
    d_size <- as.integer(floor(nn/D))
    for (d in 0:(D-1)){
        if(d != (D-1)){
            X_d <- X_matrix[(d_size*d+1):(d_size*(d+1)),,] 
            X_nod <- abind(X_matrix[0:(d_size*d),,], X_matrix[(d_size*(d+1)+1):nn,,], along=1)           
        }
        else{
            X_d <- X_matrix[(d_size*d+1):nn,,] 
            X_nod <- X_matrix[0:(d_size*d),,]           
        }
        X_nod <- convert2list(X_nod,num)
        print(class(X_nod))
        print(class(X_nod[[1]]))
        fit <- JGL(Y=X_nod, penalty=method, lambda1=lam1, lambda2=lam2, return.whole.theta=TRUE) 
        for(k in 1:num){
            sample_cov <- cov(X_d[,,k],X_d[,,k])
            tr_sum <- sum(diag(sample_cov %*% fit$theta[[i]]))
            log_det <- determinant(fit$theta[[i]], logarithm = TRUE)$modulus[1][1]
            CV_acc <- CV_acc + (tr_sum-log_det)
        }
    }
    return(CV_acc)
}
