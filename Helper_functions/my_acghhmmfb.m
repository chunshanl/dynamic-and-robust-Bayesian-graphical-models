function states = my_acghhmmfb(Pi, A, B)
%ACGHHMMFB Generate states for a Hidden Markov Model (HMM) ACGHHMMDEMO
%
%   This calculates the alpha from the HMM defined by Pi, A,and B, and
%   sample the states from alpha backwards. In Viterbi, it finds the state
%   for max(alpha) in stead.
% 
%   A(N,N) is the transition probability matrix froma state i to state j
%   B(N,T) is the emission probability matrix (B(i,k) is the probability to
%   obtain the kth log2ratio when at the state i)
%   Pi(N,1) is the distribution probability for the initial state (Pi(i) is
%   the prbability to start in the state i).
% 
%   The function return states

%   References: 
%   [1] Biological Sequence Analysis, Durbin, Eddy, Krogh, and
%       Mitchison, Cambridge University Press, 1998.
%   [2] S. Scott, "Bayesian Methods for Hidden Markov Models: Recursive
%       Computing in the 21st Century", JASA 97 (2002), pp. 337-351.


%   Copyright 2007 The MathWorks, Inc.


[numStates, T] = size(B);
states = zeros(T, 1);
alpha = zeros(numStates,T);

%------------------------------------------------------%
% Define the scaled variable for alpha
v = ones(T, 1);

% Forward
% alpha is the forward probability matrix (alpha(i,t) is the probability of
% getting to state i at time t) 
alpha(:, 1) = Pi.*B(:, 1);
v(1) = sum(alpha(:,1));
alpha(:, 1) = alpha(:,1)./v(1); % Scaled

% Scale matrix forward probability
for t = 2:T
   alpha(:, t) = A'*alpha(:, t-1) .* B(:, t);
   v(t) = sum(alpha(:,t));
   alpha(:,t) = alpha(:, t)./v(t) + realmin;
end

% Backward
states(T) = my_categorical_acghhmmsample('categorical', alpha(:, T));

for t = T-1:-1:1
    ti = (alpha(:,t) * B(:, t+1)') .* A;
    states(t) = my_categorical_acghhmmsample('categorical', ti(:, states(t+1)));
end

end % End of function

