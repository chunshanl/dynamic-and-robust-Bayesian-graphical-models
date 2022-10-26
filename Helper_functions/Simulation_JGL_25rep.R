######################
### Run fused and group lasso on three kinds of simulated data
######################
# Run Generate_and_plot_data.m and generate synthetic data
# Parallel computing on 19 cores

## Install functions 
# install.packages("R.matlab")
# install.packages("JGL")
# library(JGL)
# library(R.matlab)
# setwd("C:/1_E/2_HMM/MTHMM")
# Source functions
# source("./jgl_helper_functions.R")
# source("./get_perf_jgl.r")

## Load functions on server
# install.packages("alphashape3d", lib = 'try', dependencies = TRUE, INSTALL_opts = c('--no-lock'), repos = "http://cran.us.r-project.org")
# install.packages("rgl", lib = 'try', dependencies = TRUE, INSTALL_opts = c('--no-lock'), repos = "http://cran.us.r-project.org")
# install.packages("igraph", lib = 'try', dependencies = TRUE, INSTALL_opts = c('--no-lock'), repos = "http://cran.us.r-project.org")
# install.packages("igraph", lib = 'try', repos = "http://cran.us.r-project.org")
# install.packages("JGL", lib = 'try', dependencies = TRUE, INSTALL_opts = c('--no-lock'), repos = "http://cran.us.r-project.org")
# install.packages("JGL", lib = 'try', repos = "http://cran.us.r-project.org")
# install.packages("R.matlab", lib = 'try', dependencies = TRUE, INSTALL_opts = c('--no-lock'), repos = "http://cran.us.r-project.org")
# install.packages("foreach", lib = 'try', repos = "http://cran.us.r-project.org")
# install.packages("doSNOW", lib = 'try', repos = "http://cran.us.r-project.org")
# install.packages("iterators", lib = 'try', repos = "http://cran.us.r-project.org")
# install.packages("iterators", lib = 'try', repos = "http://cran.us.r-project.org")


################ Set up parallel computing
library(iterators, lib = 'try')
library(foreach, lib = 'try')
library(snow, lib = 'try')
library(doSNOW, lib = 'try')

# Source functions
# source("jgl_helper_functions.R")
# source("get_perf_jgl.R")
# set work directory to home
source("Helper_functions/jgl_helper_functions.R")
source("Helper_functions/get_perf_jgl.R")

ncore = 19;


cl<-makeCluster(ncore)
registerDoSNOW(cl) 
foreach(coreind= 1:ncore)  %dopar%  
# for (coreind in 1:ncore)
{

  library(igraph)
  library(JGL)
  library(R.matlab)

  ################### Set up tuning parameter grid
  ## Lambda1 options (tuning parameter for graphical lasso penalty)
  lambda1_opts <- seq(0, 0.2, by = 0.005)
  ## Lambda2 options (tuning parameter for fused or group lasso penalty)
  lambda2_opts <- c(0, 0.00005, 0.0001, seq(0.0002, 0.0096, by = 0.0004))
  ## Range of parameter values for fused graphical lasso to use when
  # computing the auc.
  lambda2s <- c(0.1, 0.15, 0.2, 0.225, 0.25, 0.35, 0.5) 
  # Note that these options were used only to obtain AUC estimates
  
  
  ################### Set up replications and run
  data_method_list = c('cmthmm', '5spikeshmm', '10spikeshmm')
  data_method_vec = rep(data_method_list, each = 25)
  rep_vec = rep(1:25, length(data_method_list))
  
  nrun = ceiling(length(data_method_list) * 25 / ncore)
  nrun_step = nrun
  if (coreind == ncore) {nrun = length(data_method_list) * 25 - (ncore-1)*nrun}

  for (jj in 1:nrun) {
    
    data_method = data_method_vec[(coreind-1)*nrun_step + jj]
    rep = rep_vec[(coreind-1)*nrun_step + jj]
    
    
    print(paste('Data: ', data_method))
    print(paste('Current rep =', rep))
    
    set.seed(10234 + rep)
    
    ################ Load true data and parameters
    cur_filename <- paste('Synthetic_data/Synthetic_data_', data_method,'_rep_', rep,'.mat', sep = '')
    data <- readMat(cur_filename)
    data <- data$data
    names(data) = dimnames(data)[[1]]
    C_true = list()
    TC = list()
    adj_true = list()
    for (s in 1:data$S.true) {
      C_true[[s]] = as.data.frame(data$C.true[,,s])
      TC[[s]] = data$TC[data$states.true == s, ]
      adj_true[[s]] =  1 * (C_true[[s]] != 0);
    }
  
    ################ Run fused joint graphical lasso
    ##Select tunning parameters
    start_time <- Sys.time()
    lambda <- get_optimal_params(TC, "fused", lambda1_opts, lambda2_opts)
    end_time1 <- Sys.time()
    end_time1 - start_time
  
    ## Fit model
    fused_res <- JGL(Y = TC, penalty = "fused", lambda1 = lambda[1], lambda2 = lambda[2],
                     return.whole.theta = TRUE)
    end_time2 <- Sys.time()
    end_time2 - start_time
  
    ## Save results
    # save estimated precision matrices
    est_C =fused_res$theta;
    C_save = c();
    for (s in 1 : data$S.true){
      C_save = rbind(C_save, est_C[[s]])
    }
    cur_theta_filename <- paste('Simulation_results/results_sim_', data_method,'_FGL_C_rep_', rep,'.csv', sep = '')
    write.table(est_C, cur_theta_filename, quote = FALSE, row.names = FALSE,
                col.names = FALSE, sep = ",")
    # save estimated adjacency matrices
    est_adj <- make.adj.matrix(fused_res$theta, separate = TRUE)
    adj_save = c();
    for (s in 1 : data$S.true){
      adj_save = rbind(adj_save, est_adj[[s]])
    }
    cur_theta_filename <- paste('Simulation_results/results_sim_',data_method,'_FGL_adj_rep_', rep,'.csv', sep = '')
    write.table(est_adj, cur_theta_filename, quote = FALSE, row.names = FALSE,
                col.names = FALSE, sep = ",")
  
    ## Get performance
    perf_summary <- get_perf_jgl(adj_true, est_adj, C_true, fused_res$theta, "fused")
    perf_summary
    # AUC
    # Run fused graphical lasso for different values of lambda2 to get range of AUC values
    # Get AUC by increasing sparsity penalty until solution is all zeros
    auc_fgl <- rep(NA, length(lambda2s))
    roc_fgl <- list()
    Y <- TC
    for (l in 1:length(lambda2s)) {
      auc_fused_res <- get_auc_jgl(param2 = lambda2s[l], adj_true, penalty = "fused")
      auc_fgl[l] <- auc_fused_res$auc
      roc_fgl[[l]] <- rbind(auc_fused_res$tpr_roc, auc_fused_res$fpr_roc)
    }
    auc_best = max(auc_fgl)
    ind_temp = which(auc_fgl == max(auc_fgl))
    roc_best = roc_fgl[[ind_temp]]
    
    

    ## Save performances
    # Take best AUC as final result + write performance summary to disk
    # tpr, fpr, mcc, auc, fl1, lambda1, lambda2, running time
    cur_perf_filename <- paste('Simulation_analysis/perf_',data_method,'_FGL_avg_rep', rep,'.csv', sep = '')
    write.table(t(c(perf_summary[1], perf_summary[2], perf_summary[3], max(auc_fgl),
                    perf_summary[7], lambda[1], lambda[2], as.numeric(end_time2 - start_time))),
                cur_perf_filename, quote = FALSE, row.names = FALSE,
                col.names = FALSE, sep = ",")
    cur_perf_filename <- paste('Simulation_analysis/perf_',data_method,'_FGL_roc_rep', rep,'.csv', sep = '')
    write.table(roc_best,
                cur_perf_filename, quote = FALSE, row.names = FALSE,
                col.names = FALSE, sep = ",")

  }

} 
stopCluster(cl)
