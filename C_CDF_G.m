%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C_CDF_G
%
%   Computes the CDF of lambda_1 of a central Wishart matrix
%   with a rank 1 mean and identity covariance. Utilizes the gamma
%   function determinant formula seen in Khatri, Kang, etc.
%
%   Inputs:      
%       M: size of Wishart matrix / Gaussian covariance parameter
%       N: degrees of freedom
%       mu1: largest eigenvalue of the Grammian of the mean matrix
%       x: domain of the CDF. Best practice is to use output from
%           ecdf of Monte Carlo data 
%
%   Returns:
%       F: CDF of lambda_1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F] = C_CDF_G(M, N, x)

 J = ones(M,1).*[1:M];
 I = J';
 
 J = J-1; %ICASSP version
 I = I-1; %ICASSP version
 
 
 F_ = @(x) det( [ gamma(N-M+I+J+1).*gammainc(x,N-M+I+J+1) ] ) / det( [gamma(N-M+I+J+1) ] ); %ICASSP version
 
% F_ = @(x) det( [ gamma(N+M-I-J+1).*gammainc(x,N+M-I-J+1) ] ) / det( [gamma(N+M-I-J+1) ] ); % Zanella version
 
 
 F = zeros(length(x),1);
 
 for k=1:length(x)
    F(k) = F_(x(k));     
 end
 