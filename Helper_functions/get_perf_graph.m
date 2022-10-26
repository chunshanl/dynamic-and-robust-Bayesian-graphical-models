function [ tpr, fpr, mcc, auc, roc_curve, fl ] = ...
    get_perf_graph( true_adj, ppi, true_C, est_C )

%% Compute performance of graph estimation averaged aross states

%% Output: 
% tpr: true positive rate of edge estimation
% fpr: false positive rate of edge estimation
% mcc: MCC, Matthews correlation coefficient of edge estimation
% auc: AUC, area under the curve of edge estimation
% roc_curve: first colomn is false positive rates and second colomn is true
% positive rates of edge estimation
% fl: Frobenius loss from comparing true and estimated precision matrix

%% Input:
% true_adj: p x p x S array of true adjacency matrices 
% ppi: p x p x S array of posterior edge inclusion probabilities
% true_C: p x p x S array of true precision matrices 
% est_C: p x p x S array of estimated precision matrices 

%%
[~, ~, K] = size(true_adj);

% 1. ---------------------------------------------------------------------
median_model = ppi > .5; % treat median model as selections

[ tp, fp, tn,fn ] = get_tp_fp_tn_fn( true_adj, median_model);

tpr = tp / (tp + fn);
fpr = fp / (fp + tn);
mcc = (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
if isnan(mcc)
    mcc=0;
end

% 2. ---------------------------------------------------------------------
opts = (0:1000) / 1000;
tpr_roc = zeros(1001, 1);
fpr_roc = zeros(1001, 1);

for i = 1:1001
    cur_threshold = opts(i);
    
    % Get estimated adjacency matrix based on ppi matrices and current threshold
    cur_adj = ppi > cur_threshold;
    
    [ tp, fp, tn,fn ] = get_tp_fp_tn_fn( true_adj, cur_adj);

    tpr_roc(i) = tp / (tp + fn);
    fpr_roc(i) = fp / (fp + tn);

end

auc = sum((fpr_roc(1:1000) - fpr_roc(2:1001)) .* ...
    (tpr_roc(2:1001) + tpr_roc(1:1000)) / 2);

roc_curve = [fpr_roc, tpr_roc];

% 3. --------------------------------------------------------------------
fl = 0;

for k = 1:K
    fl = fl + ...
        norm(squeeze(true_C(:, :, k)) - squeeze(est_C(:, :, k)), 'fro')^2 / ...
        norm(squeeze(true_C(:, :, k)), 'fro')^2 / K;
end

end

