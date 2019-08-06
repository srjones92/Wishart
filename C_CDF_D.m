%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C_CDF_D
%
%   Computes the CDF of the largest eigenvalue of a complex
%   central Wishart matrix, using the D polynomial formulation
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
function [Fd] = C_CDF_D(M, N, x)
    a = N-M;
    %a = d-1;
    
    % shift the integration limit
    %x = (x-a)/sqrt(2*a);
    
    Fd = zeros(length(x),1);
    for k=1:length(x)
        Fd(k) = Fquad(x(k),a,M);
    end
    
%     if p
%        figure(fig)
%        plot((d-1)+sqrt(2*(d-1))*x,Fd);
%        title('F_d Eq. 27');  
%        ylim([0 1]);
%     end



end

%% Distribution Computation
% Fquad
%       Quadrature integration, eq 6 ICASSP 2017 paper
function [Fd] = Fquad(x,d,M)
    a = d-1;
    D = Dpolys(a,M);
    cst = exp(-epsilon(a))/sqrt(pi);
    q = sqrt(2/a);
%    ci = 1;
%    cj = 1;
    
    Psi = zeros(M,M);
    for i=0:M-1
        ci = 1/sqrt(gamma(i+1))/sqrt(2)^i/sqrt(prod(1+(1:i)/a));
        for j=i:M-1
            cj = 1/sqrt(gamma(j+1))/sqrt(2)^j/sqrt(prod(1+(1:j)/a));
            f = @(t) polyval(D{i+1},t).*polyval(D{j+1},t).*exp(a*logt(q*t));
            intg = integral(f,-sqrt(a/2),x);
            Psi(i+1,j+1) = ci*cj*cst*intg;
            %if i==j
            %    Psi(i+1,j+1) = 1 + Psi(i+1,j+1);
            %end
            Psi(j+1,i+1) = Psi(i+1,j+1);
        end
    end
    
    Fd = det(Psi);
end

%% Helper Functions

% epsilon
%       Power series def section 5, ICASSP 2017 Paper
function [e] = epsilon(a)
    e = 1/(12*a) - 1/(360*a^3) + 1/(1260*a^5) - 1/(1680*a^7);
end

% logt
%       Def section 5, ICASSP 2017 Paper
function [s] = logt(qt)
s = log1p(qt) - qt;
%     s = zeros(size(qt));
%     for i=1:length(qt)
%         s(i) = log1p(qt(i)) - qt(i);              
%     end
end

% Dpolys
%       Generates D-poly of degree M based on recursive relations
%       Section 6 ICASSP 2017
function [D] = Dpolys(a,M)
    q = sqrt(2/a);
	D = cell(1,M);
    D{1} = 1; %D0
    D{2} = [2,-q];
    
    for n=2:M
       D{n+1} =  PolyAdd(conv([2, -(2*n+1)*q],D{n}), -(2*(n)+(n)^2*q^2)*D{n-1});
       %D{n+1}(3:end) = D{n+1}(3:end) - (2*(n-1) + (n-1)^2*(2/a))*D{n-1};
    end
end

% PolyAdd
%       Since MATLAB doesn't have a "poly" type...
function [p] = PolyAdd(p1,p2)

if iscolumn(p1)
    p1 = p1.';
end
if iscolumn(p2)
    p2 = p2.';
end

p = [zeros(1,max(0,length(p2)-length(p1))),p1] +  [zeros(1,max(0,length(p1)-length(p2))),p2];

end