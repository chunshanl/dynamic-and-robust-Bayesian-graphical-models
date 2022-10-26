function [C_save, adj_save, ppi_HMM, states_save, Phi_save, ar_Phi,...
     log_mh_ratio_Phi_save, Phi_prop_save, tau_t_save, w_kt_save, alpha_save ]=...
     MCMC_Algorithm_TDTHMM_nomean(burnin, nmc, Y, ...
     S,...      % number of hidden states 
     v0, v1, pii,...        % parameters for graphical model
     nu, a_alpha, b_alpha, K_sb, ...        % parameters for t-distribution
     Sig, C, adj, Phi,...  % initialization of graphical model
     states,...     % initialization of HMM
     tau_t, eta_kt, alpha, w_kt,...       % initialization for t-distribution
     disp_result)

% Details of inputs and outputs are in the function Call_Function_TDTHMM_nomean.m

%%
% p is number of variables; T is the number of time points 
[T, p] = size(Y);

%% Parameters and stuff in the SSSL Hao Wang prior
lambda = 1;    %  hyperparameter for the diagonal element, usually set to be 1

V0 = v0 * ones(p);
V1 = v1 * ones(p);

% Initial value for tau
% tau = repmat(V1, [1, 1, S]);  % starts with full graph
tau = repmat(V0, [1, 1, S]);
tau(adj) = v1;

ind_noi_all = zeros(p-1, p);
for i = 1:p
    if i == 1
        ind_noi = (2:p)';
    elseif i == p
        ind_noi = (1:p-1)';
    else
        ind_noi = [1:i-1,i+1:p]';
    end
    ind_noi_all(:, i) = ind_noi; 
end

ind_nok_all = zeros(S-1, S);
for i = 1:S
    if i == 1
        ind_nok = (2:S)';
    elseif i == S
        ind_nok = (1:S-1)';
    else
        ind_nok = [1:i-1,i+1:S]';
    end
    ind_nok_all(:, i) = ind_nok; 
end

%% Parameters related to tau_t
Z_it = nan(p,T);
n_kt = zeros(K_sb, T);
for t = 1:T
    for k = 1:K_sb
        Z_it( tau_t(:, t)== eta_kt(k, t), t) = k;
        n_kt(k, t) = sum(tau_t(:, t)== eta_kt(k, t));
    end
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
states_save=zeros(T, nmc);
ppi_HMM=zeros(S,T);   
n = sum(store,2)';
if disp_result
    disp(n)
end

tau_t_save = zeros(p, T, nmc);
w_kt_save = zeros(K_sb, T, nmc);
alpha_save = zeros(1, nmc);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yNorm = Y';
yNorm_tau = sqrt(tau_t).*yNorm; 
Scov=zeros(p,p,S);

% Perform MCMC sampling
    
for iter = 1: burnin + nmc
    if  mod(iter, 100) == 0 && disp_result
        fprintf('iter = %d\n', iter);
    end
    
    %% Sample tau_t
    
    b_alpha_temp = 0;
    
    for t = 1:T
        C_temp = C(:,:,states(t));
        Z_temp = Z_it(:, t);
        n_temp = n_kt(:, t);
        eta_temp = eta_kt(:,t);
        tau_temp = tau_t(:,t);
        yNorm_temp = yNorm(:,t);
        yNorm_tau_temp = yNorm_tau(:, t);
        w_temp = w_kt(:, t);
        
        % Sample z_t
        for i=1:p
            
            n_temp(Z_temp(i)) = n_temp(Z_temp(i)) -1;
            
            ind_noi = ind_noi_all(:, i);
            omega_22 = C_temp(i, i);
            omega_12 = C_temp(ind_noi, i);
            mu_c = - 1/omega_22 * omega_12' * yNorm_tau_temp(ind_noi);
            sigma_c_square = 1/ omega_22;

            temp_1 = mu_c ./ sqrt(eta_temp);
            temp_2 = sqrt(sigma_c_square ./ eta_temp);
            prob_existing = w_temp .* ...
                1./temp_2 .*( exp( - ((yNorm_temp(i) - temp_1)./ temp_2).^2./ 2  ));
            if sum(prob_existing)==0
                prob_existing = ones(length(prob_existing),1);
            end
            prob = prob_existing/sum(prob_existing);
            
            temp_1 = rand(1);
            temp_2 = [0; cumsum(prob)];
            temp_2(end) = [];
            z_ti = (temp_1 > temp_2) & (temp_1 < cumsum(prob));
            
            temp = 1:K_sb; z_ti = temp(z_ti);
            Z_temp(i) = z_ti;
            n_temp(z_ti) = n_temp(z_ti)+1;
            tau_temp(i) = eta_temp(z_ti);
            % Update yNorm_tau (x in paper)
            yNorm_tau_temp(i) = yNorm_temp(i) * sqrt(tau_temp(i));

        end
        
        Z_it(:, t) = Z_temp;
        n_kt(:, t) = n_temp;
        tau_t(:,t) = tau_temp;
        yNorm_tau(:, t) = yNorm_tau_temp;
        
        % Sample w_kt
        a_temp = 1 + n_temp(1:(K_sb-1));
        temp = cumsum(n_temp);
        b_temp = alpha + p - temp(1:(K_sb-1));
        v_vec = betarnd(a_temp, b_temp);  % sample v
        v_vec(v_vec == 1) = 0.999;
        temp = [1; cumprod(1 - v_vec)];  % compute w
        w_temp = [v_vec; 1] .* temp;
        w_kt(:, t) = w_temp; 
        
        b_alpha_temp = b_alpha_temp + sum( log(1 - v_vec  ));  % for sampling alpha
        
        % Sample eta_t 
        for k = 1:K_sb
            temp = 1:p;
            ind_k_temp = temp(Z_temp == k);
            ind_nok_temp = temp(Z_temp ~= k);
            n_k_temp = n_temp(k);
            a = (nu+n_k_temp)/2;
            b = 1/2*(nu + trace(C_temp(ind_k_temp, ind_k_temp) * yNorm_temp(ind_k_temp) *  yNorm_temp(ind_k_temp)'));
            c = trace(C_temp(ind_k_temp, ind_nok_temp) * yNorm_tau_temp(ind_nok_temp) * yNorm_temp(ind_k_temp)');
            c = c/2/sqrt(b);
            if c/sqrt(a)<= -0.7
                m = (2*a - 1) / (c + sqrt(c^2 + 4*a -2));
                u = rand(1);
                x = normrnd(m, sqrt(1/2));
                while x < 0 || m^(2*a-1)*u > x^(2*a-1) * exp(-2*(m+c)*(x-m))
                    u = rand(1);
                    x = normrnd(m, sqrt(1/2));
                end
                x = x^2;
            elseif c/sqrt(a) < 0
                m = 4*a / (sqrt(c^2 + 4*a) -c)^2;
                u = rand(1);
                b_temp = 1/m;
                x = my_gamrnd(a, b_temp);
                while u > exp(a + x * (m - 1) - 2 * sqrt(x) * c - a/m)
                    u = rand(1);
                    x = my_gamrnd(a, b_temp);
                end
            elseif c/sqrt(a) < 0.7
                u = rand(1);
                x = my_gamrnd(a, 1);
                while u > exp(-2 * c * sqrt(x))
                    u = rand(1);
                    x = my_gamrnd(a, 1);
                end
            else
                m = c + sqrt(c^2 + 4 * a);
                a_temp = 2*a;
                b_temp = 1/m; 
                u = rand(1);
                x = my_gamrnd(a_temp, b_temp);
                while u > exp(-(x - m/2 + c)^2)
                    u = rand(1);
                    x = my_gamrnd(a_temp, b_temp);
                end
                x = x^2;
            end

            eta_temp(k) = x/b;
            tau_temp(ind_k_temp) = x/b;
            % Update yNorm_tau (x in paper)
            yNorm_tau_temp(ind_k_temp) = yNorm_temp(ind_k_temp) * sqrt(x/b);

        end

        eta_kt(:,t) = eta_temp;
        tau_t(:,t) = tau_temp;
        yNorm_tau(:, t) = yNorm_tau_temp;
        
    end

    %%% Sample alpha 
    a_alpha_temp = a_alpha + T * (K_sb - 1);
    b_alpha_temp = b_alpha - b_alpha_temp;
    alpha_temp = my_gamrnd(a_alpha_temp, 1/b_alpha_temp);
    alpha = alpha_temp;
    

    %%  Construct and store centered covariance matrices for graphical models component
    for s = 1:S    
        Scov(:,:,s)=yNorm_tau(:,store(s,:)==1)*yNorm_tau(:,store(s,:)==1)';
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
            IW_scale = IW_scale + diag_nuij_inv * (e_ij) * (e_ij)' * diag_nuij_inv;
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
    if (log(rand(1)) < log_mh_ratio_Phi) || (iter<10)
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
            B(s,:) = mvnpdf( sqrt(tau_t') .* Y, zeros(1,p), Sig(:,:,s));  
        end

        B = log(B);
        B = B- repmat(max(B),S,1);
        B = exp(1).^B;
        B = .98*B + .01;
        normal = sum(B,1);
        B = B./repmat(normal,S,1);
        % Sample current states based on stationary distribution, transition matrix, and state probabilities
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
    if rem(iter,100)==0  && disp_result
        display(n)
    end

   
%% Retain values for posterior sample
    if iter > burnin
        C_save(:, :, :, iter-burnin) = C;
        adj_save(:, :, :, iter-burnin) = adj;
        Phi_save(:, :, iter-burnin) = Phi(:, :); 
        log_mh_ratio_Phi_save(iter-burnin)=log_mh_ratio_Phi;
        Phi_prop_save(:,:,iter-burnin)=Phi_prop;
        tau_t_save(:,:,iter-burnin)=tau_t;
        w_kt_save(:, :, iter-burnin) = w_kt;
        alpha_save(iter - burnin) = alpha;
    end

end %END ITERATIONS


%% Normalize ppi for classification
ppi_HMM=ppi_HMM/(nmc);


