%% Application on futures return data, run MCMC with different pi
% Parallel computing with 3 cores

mkdir('Application_results')
addpath Application_results
addpath Call_functions
addpath MCMC_algorithms
addpath Helper_functions

%%
load('futures_return_data.mat')  

%% Set up parallel computing
ncore = 3;
parpool(ncore)
% load functions 
poolobj = gcp;
addAttachedFiles(poolobj,{'Call_Function_TDTHMM_nomean_application.m','MCMC_Algorithm_TDTHMM_nomean.m',...
})

parfor coreind =1:ncore
%for coreind =1:ncore
        
    disp(coreind)
    pii_list = [2,4,5];
    pii = pii_list(coreind)/(size(data.TC,2)-1);  

    %% Run MCMC
    S = 5;  
    burnin = 2000;
    nmc = 8000;   
    h = 50^2; 
    v0 = 0.02^2; 
    v1 = h * v0;
    K_sb = 7; 
    nu = 3;
    a_alpha = 1;
    b_alpha = 1;
    
    disp_result = false;
    results=[]; states_save = [];
    rng(12351)
%     rng(12351, 'Threefry')
    [results,states_save] = Call_Function_TDTHMM_nomean_application(data, burnin, nmc, S , v0, v1, pii, nu, a_alpha, b_alpha, K_sb, disp_result);
    results.states_save = states_save;

    %% Save results
    results.states_save = [];
    tau_t_save = results.tau_t_save; results.tau_t_save = [];
    parsave(results, states_save, tau_t_save, pii_list(coreind))

end

function parsave(results, states_save, tau_t_save, pii)  %%% Set up the names here
    save(['Application_results/futures_result_pii', num2str(pii)],'results')
    save(['Application_results/futures_result_pii', num2str(pii),'_tau'],'tau_t_save')
    save(['Application_results/futures_result_pii', num2str(pii),'_states'],'states_save')
end

