%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C_CCDF_H
%
%   Computes the complementary CDF of the largest eigenvalue of a complex
%   central Wishart matrix, using the Hermite polynomial formulation
%
%   Inputs:      
%       M: size of Wishart matrix / Gaussian covariance parameter
%       N: degrees of freedom
%       x: K x 1, points at which to compute CDF
%
%   Returns:
%       Gd: CCDF at x
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Main function
function [Gd] = C_CCDF_H(M, N, x)
    %a = N-M;
    %a = d-1;
    
    % shift the integration limit
    %x = (x-a)/sqrt(2*a);

    Gd = zeros(length(x),1);
    for k=1:length(x)
        Gd(k) = Glimit(x(k),M);
    end
end

% GPsi
%       1-Fpsi, Hermitian integral iterations, section 7 ICASSP 2017 paper
function [A] = GPsi(x, M)
    A = zeros(M,M);
    A(1,1) = sqrt(pi)/2*erfc(x);
    for j=1:M-1
        A(1,j+1) = hermiteH(j-1,x)*exp(-x^2);
    end
    for i=1:M-1
        for j=i:M-1
        A(i+1,j+1) = 2*i*A(i,j) +  hermiteH(i,x)*hermiteH(j-1,x)*exp(-x^2);
        end
    end
end

% Glimit
%       Accurate computation of 1-F_limit when F_limit is near 1
function [detA] = Glimit(x,M)
    A = GPsi(x,M);
    for i=0:M-1
        for j=i:M-1
        A(i+1,j+1) = A(i+1,j+1)/sqrt(pi*2^(i+j)*factorial(i)*factorial(j));
        A(j+1,i+1) = A(j+1,i+1);
        end
    end
    if norm(A) < 0.6931471805599453 %log(2)
        detA = -detm1(A);
    else 
        detA = 1-det(eye(M)-A);
    end
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


