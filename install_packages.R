ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    
    if (length(new.pkg)){
        #chooseCRANmirror()
        install.packages(new.pkg, dependencies = TRUE,repos='http://cran.us.r-project.org')
    } 
        
    sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("latex2exp", "BDgraph", "tmvtnorm", "JGL", "corrplot", "Matrix",
              "network", "intergraph", "progress", "GGMselect", "matlab", "viridis", 
              "abind", "fields", "binr")
ipak(packages)

