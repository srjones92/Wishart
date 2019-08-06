%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   G_Threshold
%
%   Computes threshold values using the Gamma function CDF
%
%   Inputs:      
%       M: size of Wishart matrix / Gaussian covariance parameter
%       N: degrees of freedom
%       pfa: set of pfa's (use logspace)
%
%   Returns:
%       T: thresholds corresponding to pfa's
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [T] = G_Threshold(M,N,pfa)

T_G = zeros(length(pfa),1);
step_init = M;
x0 = 0;
tol = min(pfa)*10^-5;
for k  = 1:length(pfa)

    pfa_hat_G = 0;
    step_size = step_init;
    x = x0;

    iter = 0;
    while abs( pfa_hat_G - pfa(k) ) > tol 
        iter = iter + 1;
        %pfa_hat_H = 1 - wishEigCDF_H( M, N, x);
        %pfa_hat_G = 1 - C_CDF_H( M, N, x);
        pfa_hat_G = 1 - C_CDF_G( M, N, x);
        %pfa_hat_G = C_CCDF_G( M, N, x);
        if abs( pfa_hat_G - pfa(k) ) < tol || iter > 100
           T_G(k) = x;
           break;
        end

        if pfa_hat_G > pfa(k)
            x = x + step_size;
        elseif pfa_hat_G < pfa(k)
            step_size = step_size / 2;
            x = x - step_size;
        end      

    end

    T_G(k) = x;

end
%d = N - M - 1;


T = T_G;
%T = d - 1 + sqrt( 2 * (d-1) ) * T_G;

if pfa(end) == 1
   T(end) = 0; 
end

end

      %  plot( pfa, T);
       % hold on;
       % leg{ (i-1)*length(M) + j } = ['log(N) = ', num2str(log10(N(i))), ', M = ', num2str(M(j))];
        