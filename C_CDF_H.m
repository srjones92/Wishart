%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C_CDF_H
%
%   Computes the CDF of the largest eigenvalue of a complex
%   central Wishart matrix, using the Hermite polynomial formulation
%   Note: this is the distribution in the limit, not exact
%
%   Inputs:      
%       M: size of Wishart matrix / Gaussian covariance parameter
%       N: degrees of freedom
%       x: K x 1, points at which to compute CDF
%
%   Returns:
%       Fd: CDF at x
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main function
function [F] = C_CDF_H( M, N, x )
    a = N-M;
    %a = d-1;
    
    % shift the integration limit
    x = (x-a)/sqrt(2*a);
    
    
    F = zeros(length(x),1);
    for k=1:length(x)
        F(k) = Flimit(x(k),M);
    end
    
%     if p
%         figure(fig)
%         plot(x,Fd)
%     end
%     
    
end


% FPsi
%       Hermitian integral iterations, section 7 ICASSP 2017 paper
function [A] = FPsi(x,M)
    A = zeros(M,M);
    A(1,1) = sqrt(pi)/2*(1+erf(x));
    for j=1:M-1
        A(1,j+1) = -hermiteH(j-1,x)*exp(-x^2);
    end
    for i=1:M-1
        for j=i:M-1
            A(i+1,j+1) = 2*i*A(i,j) - hermiteH(i,x)*hermiteH(j-1,x)*exp(-x^2);
        end
    end
end

% Flimit
%       EQ 5, ICASSP 2017 paper
function [detW] = Flimit(x,M)
    W = FPsi(x,M);
    for i=0:M-1
        for j=i:M-1 
            W(i+1,j+1) = W(i+1,j+1)/sqrt(pi*2^(i+j)*factorial(i)*factorial(j));   
            W(j+1,i+1) = W(i+1,j+1);
        end
    end
    detW = det(W);
end