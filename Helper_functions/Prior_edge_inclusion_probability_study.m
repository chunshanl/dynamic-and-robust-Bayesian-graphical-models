%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute prior edge inclusion probabilities for different pi and v0 
%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath Call_functions
addpath MCMC_algorithms
addpath Helper_functions

%% Set up parameter grid
p = 20;
S = 3;
% x axis: pii
pii_list = [2,3,4,6,8]/(p-1);
% different lines:
v0_list = [0.02^2, 0.05^2, 0.1^2];
% fix h
h_list = 50^2;

[pii_vec, v0_vec, h_vec] = ...
    ndgrid(pii_list, v0_list, h_list);

burnin=1000;
nmc=9000;

%% Run MCMC: sample precision matrices and edges from the linked continuous spike-and-slab prior
% MCMC results are provided in the Sensitivity_results folder

% mkdir Sensitivity_results
% for ind = 1:length(pii_vec(:))
%     pii = pii_vec(ind);
%     v0 = v0_vec(ind);
%     h = h_vec(ind);
%     v1 = h * v0;
%     
%     rng(12345, 'twister')
%     results = Call_Function_LPM_sample_from_prior(burnin, nmc,  S , v0, v1, pii, p);
%     save(['Sensitivity_results/sample_from_prior_results_',num2str(ind)],'results')
%     clearvars results
% end

%% Load MCMC results and plot edge inclusion probabilities
% pii_prior = zeros(length(v0_list), length(pii_list));
% for ind = 1:length(pii_vec(:))
%     load(['Sensitivity_results/sample_from_prior_results_',num2str(ind)])
%     temp1 = ( p * (p-1) / 2 * S * nmc - sum(results.adj_save(:) == 0) / 2 ) /...
%          ( p * (p-1) / 2 * S * nmc);
%     temp2 = 1:length(v0_list);
%     rownum = temp2(v0_list == results.v0);
%     temp3 = 1:length(pii_list);
%     colnum = temp3(pii_list == results.pii);
%     pii_prior(rownum, colnum) = temp1;    
% end
% save('Sensitivity_results/pii_prior', 'pii_prior') % provided 

%%
load('Sensitivity_results/pii_prior')
linestyle = ["ro-"; "g+-"; "b*-"];
figure();
hold on;
for i = 1:length(v0_list)
    plot(pii_list, pii_prior(i,:), linestyle(i));
end
addline = refline(1, 0);
addline.Color = 'black';
set(gca,'FontSize', 13)
legend({'v_0 = 0.02','v_0 = 0.05','v_0 = 0.1', 'y = x'}, 'Location','northwest');
xlabel('$\pi$','interpreter','latex')
ylabel('$p(g_{ij}) = 1$', 'interpreter', 'latex')
hold off;
title('Prior edge inclusion probabilities')
saveas(gcf, 'Sensitivity_results/prior_edge_inclusion_probability.png')





