%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   NC_CDF_MC
%
%   Empirical estimate of the CDF of the largest eigenvalue of a complex
%   non-central Wishart matrix through Monte Carlo trials. This is
%   performed using direct Monte Carlo trials.
%
%   Inputs:      
%       M: size of Wishart matrix / Gaussian covariance parameter
%       N: degrees of freedom
%       S: N x M signal in the mean
%       nTrials: number of pseudo random trials to run
%
%   Returns:
%       F: CDF at x
%       x: domain of CDF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [F,x] = NC_CDF_MC(M, N, S, nTrials)

lambda1 = zeros(nTrials,1);
%X = zeros(N,M);

for k=1:nTrials
    X = S + 1/sqrt(2)*(randn(N,M) + 1i*randn(N,M));
    lambda1(k) = max(svd(X)).^2;    
end

[F,x] = ecdf(lambda1);



end