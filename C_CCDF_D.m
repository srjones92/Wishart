%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C_CCDF_D
%
%   Computes the complementary CDF of the largest eigenvalue of a complex
%   central Wishart matrix, using the D polynomial formulation
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
function [Gd] = C_CCDF_D( M, N, x )
    a = N-M;
    %a = d-1;
    Gd = zeros(length(x),1);
    
    % shift the integration limit
    x = (x-a)/sqrt(2*a);
    
    for k=1:length(x)
        Gd(k) = Gquad(x(k),a,M);
    end
    
%     if p
%        figure(fig)
%        semilogy(x,Gd);
%     end



end

% Gquad
%       Quadrature integration
function [Gd] = Gquad(x, d, m)    
    a = d-1;
    D = Dpolys(a,m);
    cst = exp(-epsilon(a))/sqrt(pi);
    q = sqrt(2/a);
    ci = 1;
    cj = 1;
    Psi = zeros(m,m);


    for i=0:m-1
        ci=1/sqrt(gamma(i+1))/sqrt(2)^i/sqrt(prod((1-(0:i)/a)));
        for j=i:m-1
            cj=1/sqrt(gamma(j+1))/sqrt(2)^j/sqrt(prod((1-(0:j)/a)));
            f = @(t) polyval(D{i+1},t).*polyval(D{j+1},t).*exp(a*logt(q*t));
            intg = integral(f, x, 20);       
            Psi(i+1,j+1) = ci*cj*cst*intg;

            Psi(j+1,i+1) = Psi(i+1,j+1);
        end
    end
    if norm(Psi) < 0.6931471805599453 %log(2)
        Gd =  -detm1(Psi);
    else 
        Gd = 1-det(eye(m)-Psi);
    end 
end


%% Helper Functions


% epsilon
%       Power series def section 5, ICASSP 2017 Paper
function [e] = epsilon(a)
    e = 1/(12*a) - 1/(360*a^3) + 1/(1260*a^5) - 1/(1680*a^7);
end

% logt
%       Def section 5, ICASSP 2017 Paper (as the phi function)
function [s] = logt(qt)
    s = log1p(qt) - qt;
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
