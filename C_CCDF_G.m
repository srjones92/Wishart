%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C_CCDF_G
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

function [F] = C_CCDF_G(M, N, x)

 J = ones(M,1).*[1:M];
 I = J';
 
 J = J-1; %ICASSP version
 I = I-1; %ICASSP version
 
 %C = det( gamma(N-M+I+J+1) );
 
 
 %F_ = @(x) gamma(N-M+I+J+1).*gammainc(x,N-M+I+J+1,'lower');
 %F_ = @(x)  gamma(N-M+I+J+1).*gammainc(x,N-M+I+J+1, 'lower') ;
 F_ = @(x) 1 - det( eye(M) - gamma(N-M+I+J+1).*gammainc(x,N-M+I+J+1,'lower')) / det( gamma(N-M+I+J+1) );
 %F_ = @(x) det( gamma(N-M+I+J+1).*gammainc(x,N-M+I+J+1,'lower')) / det( gamma(N-M+I+J+1) );
 %F_ = @(x) det( gamma(N-M+I+J+1).*gammainc(x,N-M+I+J+1,'upper')) / det( gamma(N-M+I+J+1) );
 F = zeros(size(x));
 
 for k=1:length(x)   
     
%     A = 1/C^(1/M)*F_(x(k));
%  
%     if norm(A) < 0.6931471805599453 %log(2)
%         F(k) = -detm1(A);
%     else 
%         F(k) = 1-det(eye(M)-A);
%     end
      F(k) = F_(x(k));     
 end
 
%  C = det( gamma(N-M+I+J+1) );
%  
%  %F_ = @(x) gamma(N-M+I+J+1).*gammainc(x,N-M+I+J+1)  ; %ICASSP version
%  
% F_ = @(x) det( [ gamma(N+M-I-J+1).*gammainc(x,N+M-I-J+1) ] ) / det( [gamma(N+M-I-J+1) ] ); % Zanella version
%  
%  
%  F = zeros(length(x),1);
%  
%  for k=1:length(x)
%     A = 1/C^(1/M)*F_(x(k));
%  
%     if norm(A) < 0.6931471805599453 %log(2)
%         F(k) = -detm1(A);
%     else 
%         F(k) = 1-det(eye(M)-A);
%     end
%  
%  end
%  
end



% detm1
%       Accurate computation of det(I-A)-1 for small A
function [detA] = detm1(A)
    tA = trace(A);
    W = A;
    pr = tA;      % previous term
    res = tA;    % result
    j = 2;
    while(abs(pr/tA) > 1.0e-16)
        W = W*A;
        pr = trace(W)/j;
        res = res + pr;
        j=j+1;
    end
    detA = expm1(-res);
end
