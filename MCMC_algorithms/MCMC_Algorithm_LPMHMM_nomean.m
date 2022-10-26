function [C_save, adj_save, ppi_HMM, states_save, Phi_save, ar_Phi,...
    log_mh_ratio_Phi_save, Phi_prop_save]...
    = MCMC_Algorithm_LPMHMM_nomean(burnin, nmc, Y, ...
    S,...    % number of hidden states 
    v0, v1, pii,...   % parameters for graphical model
    Sig, C, adj, Phi,...   % initialization of graphical model
    states,...   % initialization of HMM
    disp_result)

% Details of inputs and outputs are in the function Call_Function_LPMHMM_nomean.m

%%
% p is number of variables; T is the number of time points 
[T, p] = size(Y);

%% Parameters and stuff in the SSSL Hao Wang prior
lambda = 1;    %  hyperparameter for the diagonal element Usually set to be 1

V0 = v0 * ones(p);
V1 = v1 * ones(p);

% Initial value for tau
% tau = repmat(V1, [1, 1, S]);  % starts with full graph
tau = repmat(V0, [1, 1, S]);
tau(adj) = v1;

ind_noi_all = zeros(p-1, p);
for i = 1:p
    if i == 1
        ind_noi = [2:p]';
    elseif i == p
        ind_noi = [1:p-1]';
    else
        ind_noi = [1:i-1,i+1:p]';
    end
    ind_noi_all(:, i) = ind_noi; 
end

ind_nok_all = zeros(S-1, S);
for i = 1:S
    if i == 1
        ind_nok = [2:S]';
    elseif i == S
        ind_nok = [1:S-1]';
    else
        ind_nok = [1:i-1,i+1:S]';
    end
    ind_nok_all(:, i) = ind_nok; 
end

%% HIDDEN MARCOV MODEL Initialization
store = zeros(S,T);
for t = 1:T
    store(states(t),t)=1;
end

% Initialize the matrix with the state counts for HMM
st_counts=zeros(S,S);
% Initialize the transition matrix for states
ainit = ones(S, S);
for s = 1:S
   ainit(s,s)=ainit(s,s)+sum(store(s,:));
end
% Inital values for the HMM Stuff
A = my_dirichlet_acghhmmsample('dirichlet', ainit, S);
PI =@(x, n) (ones(1,n)/(eye(n) -x + ones(n)))';
Pi = PI(A, S);

times_of_assign = 0;

%% Posterior Storage
C_save = zeros(p, p, S, nmc);
adj_save = C_save;

log_mh_ratio_Phi_save=zeros(nmc,1);   
Phi_prop_save=zeros(S,S,nmc);
Phi_save = zeros(S, S, nmc);
% Acceptance rate for Phi
ar_Phi = 0;

% Storage for the states
states_save=zeros(T, nmc);  % zeros(T,max(nmc,500))
ppi_HMM=zeros(S,T);   
n = sum(store,2)';
if disp_result
    disp(n)
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Scov=zeros(p,p,S);

% Perform MCMC sampling


for iter = 1: burnin + nmc
    if  mod(iter, 500) == 0 && disp_result
        fprintf('iter = %d\n', iter);
    end

    
    %%  Construct and store centered covariance matrices for graphical models component
    yNorm = Y';
    for s = 1:S    
        Scov(:,:,s)=yNorm(:,store(s,:)==1)*yNorm(:,store(s,:)==1)';
    end

    %% Sample graph and precision matrix for each group using code from Hao Wang
    for s = 1:S    % state s
	    ind_nok = ind_nok_all(:, s);

        for i = 1:p   % column i
            ind_noi = ind_noi_all(:, i);

            % Equivalent to v_12 vector in Wang (2014) paper
            tau_temp = tau(ind_noi, i, s);

            % Compute m vector and diagonal entries in D matrix
            D_diag = NaN(p-1, 1);
            m = NaN(p-1, 1);
            for cur_ind = 1:p-1    % row cur_ind, column i 
                cur_ind_noi = ind_noi(cur_ind);

                % Note that stored value of tau is in fact squared version
                % Use vector of nu_ij values across S groups as diagonal elements of matrix
                % diag_nuij = diag(sqrt(squeeze(tau(cur_ind_noi, i, :))));
                tau_vec=tau(cur_ind_noi, i, :);
                diag_nuij = diag(sqrt(tau_vec(:)));

                % Cross graph relationship matrix for current edge i.e. Theta_ij in paper
                theta_xp = diag_nuij * Phi * diag_nuij;
                theta_xp11 = theta_xp(ind_nok, ind_nok);
                theta_xp12 = theta_xp(ind_nok, s);
                prod_xp = theta_xp12' / theta_xp11;
                % m(cur_ind) = prod_xp * squeeze(C(cur_ind_noi, i, ind_nok));
                C_vec=C(cur_ind_noi, i, ind_nok);
                m(cur_ind) = prod_xp * C_vec(:);
                D_diag(cur_ind) = tau_temp(cur_ind) - prod_xp * theta_xp12;
            end 
        
            Sig11 = Sig(ind_noi, ind_noi, s);
            Sig12 = Sig(ind_noi, i, s);
        
            invC11 = Sig11 - Sig12 * Sig12' / Sig(i, i, s);
        
            % Compute core of C matrix (before inversion)
            Ci = (Scov(i, i, s) + lambda) * invC11 + diag(1 ./ D_diag);
            
            % Avg with its transpose to ensure symmetry
            Ci = (Ci + Ci') ./ 2;

            % Compute Cholesky decomposition of Ci
            Ci_chol = chol(Ci);

            % Compute a vector
            a = diag(1 ./ D_diag) * m - Scov(ind_noi, i, s);

            % Use this to compute Ca where C = Ci^-1
            mu_i = Ci_chol \ (Ci_chol' \ a);

            % Generate vector from multivariate normal distribution with mean mu_i and
            % covariance matrix Ci
            % Beta here is equivalent to u vector in Wang paper
            beta = mu_i + Ci_chol \ randn(p-1, 1);
         
            % Update off-diagonal elements for current row/col
            C(ind_noi, i, s) = beta;
            C(i, ind_noi, s) = beta;
        
            % Generate sample from Gamma distribution with parameters a_gam and b_gam
            % Gam here is equivalent to v in Wang paper
            a_gam = 0.5 * n(s) + 1;
            b_gam = (Scov(i, i, s) + lambda) * 0.5;
            gam = my_gamrnd(a_gam, 1 / b_gam);
            
            % Sampled value was after change of variables. Transform to original var
            % of diagonal entry for current row/col
            c = beta' * invC11 * beta;
            C(i, i, s) = gam + c;
        
            % Updating Covariance matrix according to one-column change of precision matrix
            invC11beta = invC11 * beta;
            Sig(ind_noi, ind_noi, s) = invC11 + invC11beta * invC11beta' / gam;
            Sig12 = -invC11beta / gam;
            Sig(ind_noi, i, s) = Sig12;
            Sig(i, ind_noi, s) = Sig12';
            Sig(i, i, s) = 1 / gam;  
        
            v0 = V0(ind_noi, i);
            v1 = V1(ind_noi, i);
          
            % Bernoulli update to adjacency matrix
            w1 = -0.5 * log(v0) - 0.5 * beta .^ 2 ./ v0 + log(1 - pii);
            w2 = -0.5 * log(v1) - 0.5 * beta .^ 2 ./ v1 + log(pii);
        
            w_max = max([w1,w2], [], 2);
        
            w = exp(w2-w_max) ./ sum(exp([w1,w2] - repmat(w_max,1,2)), 2);
        
            z = (rand(p-1, 1) < w);
        
            v = v0;
            v(z) = v1(z);
        
            tau(ind_noi, i, s) = v;
            tau(i, ind_noi, s) = v;
        
            adj(ind_noi, i, s) = z;
            adj(i, ind_noi, s) = z;  
        end
    end
    
  
%% Udpate cross-graph correlation matrix Phi

    % Compute Inverse-Wishart degrees of freedom
    IW_df = p * (p - 1) / 2;
    
    % Compute Inverse-Wishart scale matrix
    V=zeros(S,S);
    for s = 1:S
        temp=C(:, :, s) .^ 2;
        temp=temp(upperInd(p));
        V(s,s) = 1 ./ sqrt(sum(temp));
    end

    IW_scale = zeros(S, S);
    for i = 1:p-1
        for j = i+1:p
            omega_ij = squeeze(C(i, j, :));
            e_ij = V * omega_ij;
            diag_nuij_inv = diag(1 ./ sqrt(squeeze(tau(i, j, :))));
            IW_scale = IW_scale + diag_nuij_inv * e_ij * e_ij' * diag_nuij_inv;
        end
    end
    
    % Sample covariance matrix
    % This gives us updated values for both Phi and V
    Psi_prop = iwishrnd(IW_scale, IW_df);
    V_inv_prop = diag(1 ./ sqrt(diag(Psi_prop)));
    Phi_prop = V_inv_prop * Psi_prop * V_inv_prop;
    
    % Compute MH acceptance ratio on log scale
    log_mh_ratio_Phi = (S + 1) / 2 * (logdet(Phi_prop) - logdet(Phi));
    
    % Accept with probability r
    if (log(rand(1)) < log_mh_ratio_Phi) || (iter==1)
        Phi = Phi_prop;
        ar_Phi = ar_Phi + 1 / (burnin + nmc);
    end


    %% MH and Gibbs sampling step for Hidden Markov Model

    % Compute the number of transitions from state i to state j
    for s = 1:S
        for r = 1:S
            st_counts(s, r) = sum((states(1:T-1) == s) .* (states(2:T) == r));
        end
    end

    % Updating block B1
    % Generate the transition matrix from the Dirichlet distributions
    Chmm = my_dirichlet_acghhmmsample('dirichlet', st_counts+1, S);
	
    % Compute the state probabilities under stationary distribution of a
    % given transition matrix C.
    PiC = PI(Chmm, S);
	
	% Create sparse representation of current state assignment
    % ind_states = sparse(states,T,1,S,T);
	
    % Compute the accepting probability using a Metropolis-Hastings step
    beta_HMM = min([1, exp(log(PiC(states(1))) - log(Pi(states(1))))]);
	
	% Accept or reject the proposed transition matrix and stationary distribution
    if rand < beta_HMM
        A = Chmm;
        Pi = PiC;
    end
    
    % Updating block B2
    % Generate copy number states using Forward propagate, backward sampling [4].
    B = zeros(S,T);
    
	% First calculate regression component based on current beta and HRF convolved X
    % Compute probability of each element belonging to each class
	% vpa function applied to increase numerical precision
	% to inrease speed of this step adjust second argument
    for s = 1:S
        B(s,:) = mvnpdf(Y, zeros(1,p), Sig(:,:,s));  
    end
    
	B = log(B);
	B = B- repmat(max(B),S,1);
	B = exp(1).^B;
    B = .98*B + .01;
	normal = sum(B,1);
	B = B./repmat(normal,S,1);
    % Sample current states based on stationary distribution, transition matrix, and stationary distribution
    states = my_acghhmmfb(Pi, A, B);  

    %%% If a state is collapsed, add observations
    if iter < burnin
      for s = 1:S
          n(s) = sum(states == s);
      end           
      for s = 1:S
          if n(s) < 21 && times_of_assign < 21
             times_of_assign = times_of_assign+1;
             if disp_result
                disp([ 'Assign: iter = ', num2str(iter), ', times = ', num2str(times_of_assign)])
                disp(n)
             end
             [num_temp, state_temp] = max(n);
             % time index of the most prevalent state
             temp = 1:T;
             temp = temp(states == state_temp);
             % first half
             temp1 = temp(1:(floor(num_temp/2)));
             % assign to state s
             states(temp1) = s;
             % update n
             n(s) = n(s) + length(temp1);
             n(state_temp) = n(state_temp) - length(temp1);
          end          
       end
    end
        
    % Given sampled states create binary indicator matrix of state classification
    for s = 1:S
        for t = 1:T 
            if(states(t)==s)
                store(s,t)=1;
				% If burn in reached add to posterior probability for each 
                if iter > burnin
                   ppi_HMM(s,t)=ppi_HMM(s,t)+1; 
                end
            else
                store(s,t)=0;
            end
        end
    end 

	% If burn in reached then add to state posterior storage
    if iter > burnin
        states_save(:,iter-burnin) = states;
    end
    
    n = sum(store,2)';
    if rem(iter,500)==0  && disp_result
        display(n)
    end
	
    
%% Retain values for posterior sample
    if iter > burnin
        C_save(:, :, :, iter-burnin) = C;
        adj_save(:, :, :, iter-burnin) = adj;
        Phi_save(:, :, iter-burnin) = Phi(:, :); 
        log_mh_ratio_Phi_save(iter-burnin)=log_mh_ratio_Phi;
        Phi_prop_save(:,:,iter-burnin)=Phi_prop;
    end

end %END ITERATIONS


%% Normalize ppi for classification
ppi_HMM=ppi_HMM/(nmc);


