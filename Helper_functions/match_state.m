function [post_order, states_est_sort] = match_state(states_true, states_est)

%% Match posterior states with true states used to generate the data 
% according to MCC of state estimation

%% Output: 
% post_order: 1 x S array of ordering of posterior states 
% states_est_sort: 1 x T array of sorted posterior states across time
%% Input:
% states_true:  1 x T array of true states 
% states_est: 1 x T array of posterior most probable state

%%
S = length(unique(states_true));

%%% list the permutation of states
states_permutation_list = perms([1:S]);
%%% get the performance of each permutation
mcc_list = zeros(size(states_permutation_list,1), 1);
for ii = 1:size(states_permutation_list,1)
    states_est_temp = states_est;
    for jj = 1:S
        states_est_temp(states_est == states_permutation_list(ii, jj)) = jj;
    end
    [~, ~, mcc] = get_perf_state(states_true, states_est_temp);
    mcc_list(ii) = mcc;
end
%%% choose the best permutation
post_order = states_permutation_list(find(mcc_list == max(mcc_list)), :);

states_est_sort = states_est;
for s = 1:S
    states_est_sort(states_est == post_order(s)) = s;
end



