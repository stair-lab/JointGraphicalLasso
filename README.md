# Joint Gaussian Graphical Model Estimation: A Survey

This repo is the implementation of the review paper [[1]](#1)
## Test Models
1. Fused graphical lasso [[2]](#2)
2. Group graphical lasso [[2]](#2)
3. Graphical lasso [[2]](#2)
4. Doubly joint spike-and-slab graphical lasso [[3]](#3)
    
## Installation
1. Anaconda Environment package: 
```
conda env create -f environment.yml
conda activate r_env2  #activate environment
```
2. Install R packages
```
Rscript install_packages.R
```

## Run Examples

#### Jupyter notebook
Saveral examples of data generation processes as well as sample codes are in the folder `./examples/jupyter_notebook`

#### Plot ROC curve
Sample code for data generation process 1 (DGP1). The instruction for running DGP2_roc.r is the same.

```
cd examples/roc
### Generate simulated data, the result will be stored in ./data 
Rscript DGP1_roc.r DG [DATA DIMENSION]

### Select one of the refularization method FGL/GGL/GL. The result will be stored in ./results
Rscript DGP1_roc.r [ACTION: FGL/DGL/GL] [DATA DIMENSION]

###visualization
Rscript DGP1_roc_visualization.r
```


##### Other examples
Please check the structure tree below for more details.


## Structure
```bash
├── examples
│   ├── jupyter_notebook
|   |   ├── simple_example_block.ipynb
|   |   ├── simple_example_scalefree.ipynb
|   |   ├── simple_example_ssjgl.ipynb
│   │   └── simple_example.ipynb
│   │
│   ├── roc # run & visualize ROC curve
|   |   ├── DGP1_roc_visualization.r #visualization
│   |   ├── DGP1_roc.r # roc curve on scalefree network, common structures share same inverse convarince matrix (data generation process 1)
|   |   |                
|   |   ├── DGP2_roc_visualization.r #visualization
|   |   ├── DGP2_roc.r # roc curve on scalefree network, common structures have different inverse convarince matrices (data generation process 2)
|   |   |                    
|   |   ├── simple_roc_vis.r # visualization
|   |   └── simple_roc.r # roc curve on ramdom network
|   | 
|   ├── joint_demo.r # beautiful result on random network (Erdos-Renyi graph)            
│   ├── loss_graphsize_npAIC.r #fix p, vary n            
│   ├── loss_smallgraphsize.r #fix n, vary n             
│   ├── oos_scalefree.r # out-of-sample likelihood on scalefree network.              
│   ├── oos.r # out-of-sample likelihood on random network      
|   ├── scalefree_AIC.r # model selection on scalefree network using AIC, tune the trucation value                
|   ├── scalefree_BIC.r # model selection on scalefree network using BIC, tune the trucation value               
|   ├── simple_example_ar.r # example on AR network: model selction, fnr,fpr, Frobenious loss, etropy loss                      
|   └── simple_example_scalefree.r # example on scalefree network: model selction, fnr,fpr, Frobenious loss, etropy loss
|                          
├── R #source file
|   ├── admm.iters.R
|   ├── display.R
|   ├── eval.R
|   ├── gen_data.R
|   ├── gete.R
|   ├── JGL.R
|   ├── metrics.R
|   └── SSJGL.R
|   
├── environment.yml
├── install_packages.R
├── README.md
└── .gitignore
```

## References

<a id="1">[1]</a> 
Tsai, K., Koyejo, O., & Kolar, M. (2022). 
Joint Gaussian graphical model estimation: A survey. 
Wiley Interdisciplinary Reviews: Computational Statistics, e1582.

<a id="2">[2]</a> 
Danaher, P., Wang, P., & Witten, D. M. (2014). 
The joint graphical lasso for inverse covariance estimation across multiple classes. 
Journal of the Royal Statistical Society. Series B, Statistical methodology, 76(2), 373.

<a id="3">[3]</a>
Zehang Richard Li, Tyler H. McCormick, and Samuel J. Clark. 
["Bayesian joint spike-and-slab graphical lasso"](https://github.com/richardli/SSJGL).
 International Conference on Machine Learning, 2019.

