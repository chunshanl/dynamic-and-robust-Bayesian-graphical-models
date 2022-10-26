function [tpr, fpr, mcc] = get_perf_state(states_true, states_est)

%% Compute performance of state estimation averaged aross states
%% Output: 
% tpr: true positive rate
% fpr: false positive rate
% mcc: Matthews correlation coefficient
%% Input:
% states_true:  1 x T array of true states 
% states_est: 1 x T array of posterior most probable state

%%
S = length(unique(states_est));

tp = 0;
tn = 0;
fp = 0;
fn = 0;

for s = 1:S
    temp_est = states_est == s;
    temp_true = states_true == s;
    tp = tp + sum(temp_est & temp_true);
    tn = tn + sum((~temp_est) & (~temp_true));
    fp = fp + sum((temp_est) & (~temp_true));
    fn = fn + sum((~temp_est) & (temp_true));
end

tpr = tp / (tp + fn);
fpr = fp / (fp + tn);
mcc = (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));

