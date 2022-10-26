# Dynamic and Robust Bayesian Graphical Models
_____________________________

Author: Chunshan Liu

Contact: chunshanl@hotmail.com

This is the Matlab and R code to reproduce results in the following paper:
>Dynamic and Robust Bayesian Graphical Models, Chunshan Liu, Daniel R. Kowal, and Marina Vannucci
______________________________

## Introduction

The repository contains the Matlab code for the proposed Dynamic and Robust Bayesian Graphical Model. The proposed model is a Bayesian graphical model for heavy-tailed time series data. It provides interpretable representations and insightful visualizations of the relationships among time series. For dynamics, the model employ state-of-the-art hidden Markov models (HMMs) and outputs hidden states with state-specific graphs. In addition, it generates parameters that measure the similarity between each pair of states. Last, the model generates the volatility estimation of each time series across time. 

An example of applying the proposed model and doing posterior analysis can be found in Simulation_study.m. The repository also contains Matlab and R code to reproduce all results in the simulation study, sensitivity study and application study.

## Files:

- Simulation_study.m
  - This script can reproduce MCMC results, plots and tables in the simulation study. It includes a script to randomly generate precision matrices and data sets. The data and results are exactly reproducible from the scripts provided since the random number generator seeds have been fixed.
  - This script contains instructions to run the three proposed Bayesian graphical models: the dynamic Gaussian graphical model, the dynamic classical-t graphical model and the dynamic Dirichlet-t graphical model. Search for section "Run MCMC: an example".
  - This script demonstrates functions to plot MCMC results, generate model performances, and perform MCMC diagnostics as discribed in the paper.

- Sensitivity_study.m
  - This script can reproduce results in the sensitivity study. 

- Application_gesture_return.m
  - This script can reproduce results in the application study on the hand gesture data. 

- Call_functions
  - This folder contains functions to call the MCMC algorithms of the proposed Bayesian models in this paper. These functions generate initial values, and pass data and hyperparameters into the MCMC algorithms.

- MCMC_algorithms
  - This folder contains the MCMC algorithms of proposed dynamic Bayesian Graphical models.

- Helper functions
  - This folder contains helper functions for the MCMC algorithms and for regenerating the results presented in the paper. They are called in the script and don't need to be run independently.

- Application_data
  - This folder contains the Gesture Phase Segmentation Data Set available on the [UCI Machine Learning Repository](https://archive.ics.uci.edu/ml/index.php).

## Acknowledgements

The code provided here either includes or calls code associated with the following publications:

>Warnick, R., Guindani, M., Erhardt, E., Allen, E., Calhoun, V., and Vannucci, M. (2018).A bayesian approach for estimating dynamic functional network connectivity in fmri data.Journal of the American Statistical Association, 113(521):134–151

>Peterson,  C.  B.,  Osborne,  N.,  Stingo,  F.  C.,  Bourgeat,  P.,  Doecke,  J.  D.,  and  Vannucci,M.  (2020).   Bayesian  modeling  of  multiple  structural  connectivity  networks  during  theprogression of alzheimer’s disease. Biometrics.

> Wang H. (2015). Scaling it up: Stochastic search structure learning in graphical models. Bayesian Analysis, 10(2), 351-377.