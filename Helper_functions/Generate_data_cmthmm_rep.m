%% Generate classical-t data for 25 iterations
%% Parameters
p = 20; % number of variables
S_true = 3; % number of hidden states
T = S_true*250;  % length of the time series
rng(12345, 'twister');

nblock=15;  % divide the time series equally into nblock number of blocks
n_perturb_delete = [10,10];  % number of edges to delete from s-1 to s
n_perturb_add = [10,10];  % number of edges to add from s-1 to s
C_true=zeros(p,p,S_true);
nu_true=3; % true value of nu

%% Generate precision matrices
% AR(2) graph in state 1
C_true(:,:,1) = toeplitz([1, 0.5, 0.4, zeros(1, p - 3)]);
C_temp=C_true(:,:,1);

for k=2:S_true
    
    % Locations of all nonzero entries above the diagonal
    [rpos, cpos] = find(triu(C_temp) - eye(p));

    % Locations of all zero entries above the diagonal
    [rzero, czero] = find(triu(C_temp == 0));

    % Sample locations for zeros and nonzero entries without replacement
    % to perturb
    pos_inds = randsample(size(rpos, 1), n_perturb_delete(k-1));
    zero_inds = randsample(size(rzero, 1), n_perturb_add(k-1));

    % Add new edges
    for j = 1:n_perturb_add(k-1)
        % Generate random sign for new entry
        sign = 1;
        if binornd(1, .5) == 0
            sign = -1;
        end

        % Generate nonzero value for precision matrix entry
        nonzero_val = unifrnd(.4, .6) * sign;
        C_temp(rzero(zero_inds(j)), czero(zero_inds(j))) = nonzero_val;
        C_temp(czero(zero_inds(j)), rzero(zero_inds(j))) = nonzero_val;
    end   
    % Delete existing edges
    for j = 1:n_perturb_delete(k-1)
        % Assign a current nonzero value to be 0
        C_temp(rpos(pos_inds(j)), cpos(pos_inds(j))) = 0;
        C_temp(cpos(pos_inds(j)), rpos(pos_inds(j))) = 0;
    end

    % Adjust revised matrices to ensure positive definiteness
    C_true(:,:,k) = fix_matrix(C_temp, 1);
    
    % Verify that results are positive definite
    all(eig(C_true(:,:,k)) > 0);
end

%% Compute Phi
Phi = zeros(S_true, S_true);
for i = 1:(S_true - 1)
    for j = (i+1):S_true
        a = C_true(:,:,i);
        a = a(upperInd(p));
        b = C_true(:,:,j);
        b = b(upperInd(p));
        Phi(i,j) = corr(a,b);
    end
end
data.Phi_true=Phi;

%% Generate hidden states 
rng(3, 'twister')
states_true = zeros(1, T);
s = randsample(S_true, nblock, true);  % replacement=true; randomly assign states to nblock number of blocks
blocklength = floor(T/nblock);
for i = 1:nblock
    a = ((i - 1)*blocklength + 1):i * blocklength;
    states_true(a) = s(i);
end
states_true(i * blocklength:T) = s(end);

%% Generate Data for 25 replications
nrep =25;
for rep = 1:nrep    

    rng(12345 + rep, 'twister');
    data = [];
    
    %% Generate MVN data
    TC = zeros(T,p);
    for k = 1:S_true
        TC(states_true == k, :) = mvnrnd(zeros(p, 1), inv(C_true(:,:,k)), sum(states_true == k));
    end
    
   %% Transform MVN data to Classical t-distribution data
    tau_t_true = gamrnd(nu_true/2, 2/nu_true, [T, 1]);
    TC_new = TC ./ sqrt(tau_t_true);
    TC = TC_new;

    %% Save data
    data.TC=TC;
    data.S_true=S_true;
    data.C_true=C_true;
    data.states_true=states_true;
    data.Phi_true=Phi;
    data.nu_true=nu_true;
    data.tau_t_true=tau_t_true;
    data.method = 'cmthmm';
    save(['Synthetic_data/Synthetic_data_cmthmm_rep_', num2str(rep),'.mat'],'data')

end



 