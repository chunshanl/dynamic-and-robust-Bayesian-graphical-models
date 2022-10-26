function outval = my_dirichlet_acghhmmsample(distrtype, varargin)
%ACGHHMMSAMPLE The random sampling function for ACGHHMMDEMO
%
%   This function return randon sampling from 
%       truncated normal distributions - normal, type = 1
%       truncated gamma distributions - gamma, type = 2
%       dirichlet distribution - dirichlet, type = 3
%       categorical distribution - categorical, type = 4
% 
%   Input for the truncated normal distribution: 
%           mu, tau, lower-bound and upper-bound.
% 
%   Input for the truncated gamma distribution: 
%           alpha, beta, state.
% 
%   Input for the dirichlet distribution: 
%           a row vector, number of rows in return matrix.
%
%   Input for the categorical distribution:
%           a vector of distribution probability

%   Copyright 2007-2009 The MathWorks, Inc.


distributions = {'normal','gamma','dirichlet','categorical'};
% k = find(strncmpi(distrtype, distributions, numel(distrtype)));
% outval = [];
% switch(k)
%     case 1  % normal - compute mu from prior mu normal distribution
%         truncnormrand = @(a,b) norminv(normcdf(a) + rand*(normcdf(b)-normcdf(a)));
%         mu = varargin{1};
%         tau = varargin{2};
%         lowerbound = (varargin{3} - mu)/tau;
%         upperbound = (varargin{4} - mu)/tau;
%         
%         outval = truncnormrand(lowerbound, upperbound)*tau + mu;
%     case 2 % gamma - compute sigma from gamma distribution
%         alpha = varargin{1};
%         beta = varargin{2};
%         bound = (varargin{3})^2;
%         
%         x = sampleTruncatedGamma(alpha, 1/beta, 1/bound);
%         outval = sqrt(1/x);       
%     case 3 % dirichlet
        a = varargin{1}; % A row vector
        
        ad = gamrnd(a, 1);
        ad_rowsum = sum(ad, 2);
        outval = ad ./ repmat(ad_rowsum, 1, size(a, 2));
%     case 4 % categorical
%         p = varargin{1};
%         cdf = cumsum(p);
%         if cdf(end) <=0
%             return;
%         end
%         outval = sum(cdf < rand*cdf(end)) + 1; 
% end     
end % End of function
% 
% function x = sampleTruncatedGamma(a, b, lowbound)
% % x -random number from truncated gamma with shape and scale parameters A
% % and B and x > lowbound
% 
% % Reference: "Non-Uniform Random Variate Generation" L.Devroye (1986)
% low = lowbound/b;
% if a < 1
%     x = low + exprnd(1);
%     u = rand;
%     while (log(x)+log(u)/(1-a)) > log(a)
%        x = low + exprnd(1);
%        u = rand;
%     end
% elseif a== 1
%     x = low + exprnd(1);
% else
%     if low > a-1
%         e1 = exprnd(1);
%         e2 = exprnd(2);
%         x = low + e1/(1-(a-1)/low);
%         while (x/low-1+log(low/x)) > (e2/(a-1))
%             e1 = exprnd(1);
%             e2 = exprnd(1);
%             x = low + e1/(1-(a-1)/low);
%         end
%     else
%         x = gamrnd(a, 1);
%         while x <= low
%             x = gamrnd(a, 1);
%         end
%     end
% end
% x = x*b;
% end % End of function
