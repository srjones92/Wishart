%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   NC_CDF_G
%
%   Computes the CDF of lambda_1 of a non-central Wishart matrix
%   with a rank 1 mean and identity covariance. Note that the
%   normalizing constant is not explicity, rather the maximum value
%   is used. 
%
%   Inputs:      
%       M: size of Wishart matrix / Gaussian covariance parameter
%       N: degrees of freedom
%       mu1: largest (only non-zero) eigenvalue of the Grammian 
%           of the mean matrix
%       x: domain of the CDF. Best practice is to use output from
%           ecdf of Monte Carlo data 
%
%   Returns:
%       F: CDF of lambda_1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F] = NC_CDF_G(M, N, mu1, x)

 J = ones(M,1).*[1:M];
 I = J';
 
 % Hypothesized Correct Version - Derived from Cauchy Binet
% F_ = @(x) det( [ integral( @(t) hypergeom([],N-M+1,mu1.*t).*t.^(N-I(:,1)).*exp(-t), 0, x,'ArrayValued',true), ...
%     gamma(N+M-I(:,2:end)-J(:,2:end)+1).*gammainc(x,N+M-I(:,2:end)-J(:,2:end)+1) ] );
 
 
 
% Test versions of function

% Normalization from Zanella 2009 eq 6
 F_ = @(x) det( [ integral( @(t) hypergeom([],N-M+1,mu1.*t).*t.^(N-I(:,1)).*exp(-t), 0, x,'ArrayValued',true), ...
     gamma(N+M-I(:,2:end)-J(:,2:end)+1).*gammainc(x,N+M-I(:,2:end)-J(:,2:end)+1) ] );
% No change - multilinear map, 1/factorial(N-M) just folds into the
% constant. Why introduce it in the paper?
 
 
% Flipped-Indices
% F_ = @(x) det( [ integral( @(t) hypergeom([],N-M+1,mu1.*t).*t.^(N-I(:,1)).*exp(-t), 0, x,'ArrayValued',true), ...
%     gamma(N-M+I(:,2:end)+J(:,2:end)+1).*gammainc(x,N-M+I(:,2:end)+J(:,2:end)+1) ] );
 
 F = zeros(length(x),1);
 
for k=1:length(x)
    F(k) = F_(x(k));    
end


% implicit normalization
F = F./max(F);

end