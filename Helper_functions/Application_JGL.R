######################
### Run fused lasso on gesture data given states estimated by our proposed model
######################

library(igraph)
library(JGL)
library(R.matlab)

# Source functions
source("Helper_functions/jgl_helper_functions.R")
source("Helper_functions/get_perf_jgl.R")
  
## Load data -------------
cur_filename <- 'Application_data/a1_va3.mat'
data <- readMat(cur_filename)
p = dim(data$TC)[2]
data$TC_standardized = data$TC
for (i in 1:p){
  data$TC_standardized[,i] = data$TC[, i] / sd(data$TC[, i])
}
# load estimated states from the proposed model
states_dirichlet = read.csv('states_dirichlet.csv', header = F)
S = dim(table(states_dirichlet))
# shape data for the JGL function
TC = list()
for (s in 1:S) {
  TC[[s]] = data$TC_standardized[states_dirichlet == s, ]
}

## Run fused joint graphical lasso -------------------------
## Lambda1 options (tuning parameter for graphical lasso penalty)
lambda1_opts <- seq(0, 0.01, by = 0.0005)
#lambda1_opts <- seq(0, 0.2, by = 0.005)
## Lambda2 options (tuning parameter for fused or group lasso penalty)
lambda2_opts <- c(0, 0.00005, 0.0001, seq(0.0002, 0.0096, by = 0.0004))

set.seed(12345)
## Select tunning parameters based on AIC 
# The result is [0, 0] and the fitted model generates full graphs
# lambda <- get_optimal_params(TC, "fused", lambda1_opts, lambda2_opts)

## Fit model, the best lambda that matches our proposed model is [0.1, 0.09]
fused_res <- JGL(Y = TC, penalty = "fused", lambda1 = 0.1, lambda2 = 0.09,
                 return.whole.theta = TRUE)

## Save results ------------------------------
# save estimated precision matrices
est_C =fused_res$theta;
C_save = c();
for (s in 1 :S){
  C_save = rbind(C_save, est_C[[s]])
}
cur_theta_filename <- 'Application_results/gesture_result_FGL_C.csv'
write.table(C_save, cur_theta_filename, quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = ",")
# save estimated adjacency matrices
est_adj <- make.adj.matrix(fused_res$theta, separate = TRUE)
# for (s in 1:S){
#   print((sum(est_adj[[s]]) - p)/2)
# }
adj_save = c();
for (s in 1 : S){
  adj_save = rbind(adj_save, est_adj[[s]])
}
cur_theta_filename <- 'Application_results/gesture_result_FGL_adj.csv'
write.table(adj_save, cur_theta_filename, quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = ",")

