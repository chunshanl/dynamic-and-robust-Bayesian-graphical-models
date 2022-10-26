% Application on gesture data, run MCMC with different S
% Parallel computing with 5 cores

mkdir('Application_results')
addpath Application_results
addpath Call_functions
addpath MCMC_algorithms
addpath Helper_functions

%%
load('Application_data/a1_va3.mat')
data.TC = TC;
TC = [];

%% Set up parallel computing
ncore = 4;
parpool(ncore)
% load functions 
poolobj = gcp;
addAttachedFiles(poolobj,{'Call_Function_TDTHMM_nomean_application.m','MCMC_Algorithm_TDTHMM_nomean.m',...
})

parfor coreind =1:ncore
%for coreind =1:ncore
        
    disp(coreind)

    %% Run MCMC
    S = coreind + 1;  %%% Set up the number of states here
    burnin = 8000;
    nmc = 8000;   
    h = 50^2; 
    v0 = 0.02^2; 
    v1 = h * v0;
    pii = 3/(size(data.TC,2)-1);  
    K_sb = 6; 
    nu = 3;
    a_alpha = 1;
    b_alpha = 1;
    
    disp_result = false;
    results=[]; states_save = [];
    rng(12345, 'Threefry')
    [results,states_save] = Call_Function_TDTHMM_nomean_application(data, burnin, nmc, S , v0, v1, pii, nu, a_alpha, b_alpha, K_sb, disp_result);
    results.states_save = states_save;

    %% Save results
    results.states_save = [];
    tau_t_save = results.tau_t_save; results.tau_t_save = [];
    parsave(results, states_save, tau_t_save, S)

end

function parsave(results, states_save, tau_t_save, S)  %%% Set up the names here
    save(['Application_results/gesture_result_', num2str(S),'states'],'results')
    save(['Application_results/gesture_result_', num2str(S),'states_tau'],'tau_t_save')
    save(['Application_results/gesture_result_', num2str(S),'states_states'],'states_save')
end

