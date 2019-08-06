%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C_CDF_MC
%
%   Empirical estimate of the CDF of the largest eigenvalue of a complex
%   central Wishart matrix, through Monte Carlo trials generated using the
%   Bartlett decomposition function wishrndC
%
%   Inputs:      
%       M: size of Wishart matrix / Gaussian covariance parameter
%       N: degrees of freedom
%       nTrials: number of pseudo random trials to run
%
%   Returns:
%       F: CDF at x
%       x: domain of CDF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [F,x] = C_CDF_MC(M, N, nTrials)

lambda1 = zeros(nTrials,1);

for k=1:nTrials
    lambda1(k) = max(eig( wishrndC(eye(M),N,eye(M),1) ));
end

[F,x] = ecdf(lambda1);



end