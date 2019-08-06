%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   D_Threshold
%
%   Computes threshold values using wishEigCDF_D
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

function [T] = D_Threshold(M,N,PF)

T_D = zeros(length(PF),1);
step_init = M;
x0 = 0;
tol = min(PF)*10^-5;

max_iter = 200;

 for k  = 1:length(PF)

    pfa_hat_D = 0;
    step_size = step_init;
    x = x0;

    iter = 0;
    while abs( pfa_hat_D - PF(k) ) > tol 
        iter = iter + 1;
        
        %pfa_hat_D = C_CCDF_D( M, N, x);
        pfa_hat_D =  1- C_CDF_D( M, N, x);
        if abs( pfa_hat_D - PF(k) ) < tol || iter > max_iter
           T_D(k) = x;
           break;
        end

        if pfa_hat_D > PF(k)
            x = x + step_size;
        elseif pfa_hat_D < PF(k)
            step_size = step_size / 2;
            x = x - step_size;
        end      

    end

    T_D(k) = x;

 end
a = N - M ;

T = T_D;
T = a + sqrt( 2 * (a) ) * T_D;
if PF(end) == 1
   T(end) = 0; 
end

end
