%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   H_Threshold
%
%   Computes threshold values using wishEigCDF_H
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


function [T] = H_Threshold(M,N,pfa)

T_H = zeros(length(pfa),1);
step_init = M;
x0 = 0;
tol = min(pfa)*10^-3;
for k  = 1:length(pfa)

    pfa_hat_H = 0;
    step_size = step_init;
    x = x0;

    iter = 0;
    while abs( pfa_hat_H - pfa(k) ) > tol 
        iter = iter + 1;
        %pfa_hat_H = 1 - wishEigCDF_H( M, N, x);
        %pfa_hat_H = 1 - C_CDF_H( M, N, x);
        pfa_hat_H = C_CCDF_H( M, N, x);
        if abs( pfa_hat_H - pfa(k) ) < tol || iter > 100
           T_H(k) = x;
           break;
        end

        if pfa_hat_H > pfa(k)
            x = x + step_size;
        elseif pfa_hat_H < pfa(k)
            step_size = step_size / 2;
            x = x - step_size;
        end      

    end

    T_H(k) = x;

end
d = N - M - 1;

T = d - 1 + sqrt( 2 * (d-1) ) * T_H;

if pfa(end) == 1
   T(end) = 0; 
end

end

      %  plot( pfa, T);
       % hold on;
       % leg{ (i-1)*length(M) + j } = ['log(N) = ', num2str(log10(N(i))), ', M = ', num2str(M(j))];
        