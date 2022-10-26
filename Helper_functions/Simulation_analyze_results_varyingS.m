%% Compute WAIC for each replication
S = 2;
data_method = "10spikeshmm";
model_method = "TDTHMM";
waic_vec = zeros(1, 25);

for rep = 1:25

    data = [];
    data = load(strcat('Synthetic_data/Synthetic_data_', data_method, '_rep_', num2str(rep)));
    data = data.data;
    [T, p] = size(data.TC);
    % when S is not 3
    load(strcat('Simulation_results_varyingS/sim_results_',data_method,'_',model_method, '_', num2str(S),'states_rep_', num2str(rep)))
    load(strcat('Simulation_results_varyingS/sim_results_',data_method,'_',model_method, '_', num2str(S), 'states_states_rep_', num2str(rep)))
    % when S is 3
    %load(strcat('Simulation_results/sim_results_',data_method,'_',model_method,'_rep_', num2str(rep)))
    %load(strcat('Simulation_results/sim_results_',data_method,'_',model_method,'_states_rep_', num2str(rep)))

    results.states_save = states_save;

    %% Compute WAIC
    TC_scaled = data.TC;
    if isfield(results, 'mean_input')
        TC_scaled = TC_scaled - results.mean_input;
    end
    if isfield(results, 'std_input')
        TC_scaled = TC_scaled./results.std_input;
    end
    log_likelihoods = zeros(results.nmc, T);
    if strcmp(results.method, 'TDTHMM')
        for t = 1:T
            for ii = 1:results.nmc
                C_temp = results.C_save(:,:,results.states_save(t, ii), ii);
                %adj_temp = ~results.adj_save(:,:,results.states_save(t, ii),ii);
                %C_temp(adj_temp) = 0;
                C_temp = diag(sqrt(results.tau_t_save(:,t,ii))) * C_temp * diag(sqrt(results.tau_t_save(:,t,ii)));
                temp = -p/2 * log(2*pi) + 1/2 * log(det(C_temp)) - ...
                        1/2 * TC_scaled(t,:) * C_temp * TC_scaled(t,:)';
                log_likelihoods(ii, t) = temp;
    %             if ~isreal(temp)
    %                 log_likelihoods(ii, t) = 0;
    %             end
            end
        end
    elseif  strcmp(results.method, 'CMTHMM')
        for t = 1:T
            for ii = 1:results.nmc
                C_temp = results.C_save(:,:,results.states_save(t, ii), ii);
                % adj_temp = ~results.adj_save(:,:,results.states_save(t, ii),ii);
                % C_temp(adj_temp) = 0;
                C_temp = results.tau_t_save(t,ii)* C_temp;
                temp = -p/2 * log(2*pi) + 1/2 * log(det(C_temp)) - ...
                        1/2 * TC_scaled(t,:) * C_temp * TC_scaled(t,:)';
                log_likelihoods(ii, t) = temp;
    %             if ~isreal(temp)
    %                 log_likelihoods(ii, t) = 0;
    %             end
            end
        end
    end
    
    lppd = sum(log(mean(exp(log_likelihoods))));
    pwaic2 = sum(var(log_likelihoods));
    waic = - 2 * lppd + 2 * pwaic2;
    % disp(waic)

    %% Save
    waic_vec(rep) = waic;

end