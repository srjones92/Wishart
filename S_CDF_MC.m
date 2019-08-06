%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   S_CDF_MC
%
%   Empirical estimate of the CDF of the largest eigenvalue of a complex
%   central Wishart matrix with arbitrary covariance matrix, through Monte 
%   Carlo trials generated using the Bartlett decomposition function 
%   wishrndC
%
%   Inputs:      
%       M: size of Wishart matrix / Gaussian covariance parameter
%       N: degrees of freedom
%       Sigma: Covariance matrix
%       nTrials: number of pseudo random trials to run
%
%   Returns:
%       F: CDF at x
%       x: domain of CDF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [F,x] = S_CDF_MC(M, N, Sigma, nTrials)

lambda1 = zeros(nTrials,1);
D = chol(Sigma);

for k=1:nTrials
    lambda1(k) = max(eig( wishrndC(Sigma,N,D,1) ));
end

[F,x] = ecdf(lambda1);



end