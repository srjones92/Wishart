%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   NC_CDF_D
%
%   Computes the CDF of lambda_1 of a non-central Wishart matrix
%   in the MxM case, using Laguerre polynomial inner products from the
%   central case, albeit with the hypergeometric terms appearing in the
%   first column. 
%
%   Inputs:      
%       N: degrees of freedom
%       M: size of matrix
%       mu1: largest (only non-zero) eigenvalue of the Grammian 
%           of the mean matrix
%       x: domain of the CDF. Best practice is to use output from
%           ecdf of Monte Carlo data 
%
%   Returns:
%       F: CDF of lambda_1
%       Psi_: Matrices corresponding to each element of x
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F, Psi_] = NC_CDF_D_Backup( M, N, mu1, x ) %nc_wishEig_CDF_Laguerre( N, M, mu1, x ) 
    % 
    
    F01 = @( t ) hypergeom( [], N-M+1, mu1*t );
    
    F = zeros(length(x),1);
  
    a = N - M;

    
    % Generate the Laguerre coefficients
    L = LaguerreCoef(M,a);
    
    % Generate the D coefficients
    D = Dpolys( a, M );

    %upper limit shift
    %y = x/sqrt(2*a) - sqrt(a/2);
    y = (x-a)/sqrt(2*a);
    
    % compute normalization coef
    Flim = Normalization(N, M, mu1, y, D);
        
    
    
    for ind=1:length(x)
    
        Psi = zeros(M,M);
        
        for i = 1:M 
         
            %C_i1 = cij(a,i,1) * a^(-i/2) * a^(-1/2)*sqrt(factorial(i)*factorial(1))/sqrt(pi);
            C_i1 = cij(a,i,1) / sqrt( pi * 2^(i+1) * factorial(i) ); %j=1
            
            for l=i:M      
                H = 0;
                for k=2:M
                    H = H + integral( @(t) exp( phi(a,t) ) .* D{M}(k) .* D{M+1-i}(l-i+1) .* t.^(2*M-l-k), ...
                        -sqrt(a/2), y(ind) );
                end
                
                Psi(i,1) = Psi(i,1) + C_i1 * ( H + integral( @(t) exp( phi(a,t) ) .* D{M}(1) .* D{M+1-i}(l-i+1) .* ...
                    F01( a + t*sqrt(2*a) ) .* t.^(M-l), -sqrt(a/2), y(ind) ) );

            end
                    

            for j=2:M 
               Psi(i,j) = cij(a,i,j) / sqrt( pi * 2^(i+j) * factorial(i) * factorial(j) ) * ...
                    integral( @(t) exp( phi(a,t) ) .* polyval( D{M-i+1}, t ) ...
                    .* polyval( D{M-j+1}, t ), -sqrt(a/2), y(ind) );

            end
        end
        
%        Testing the Stephen methodology...
%         if norm(Psi) < 0.6931471805599453 %log(2)
%             F(ind) =  -detm1(Psi);
%         else 
%             F(ind) = 1-det(eye(M)-Psi);
%         end 
     
% what value of mu1 pushes cond(Psi) > 1/eps
        if cond(Psi) > 1/eps
            F(ind) = det(Psi);
        else
            F(ind) = abs(prod(eig(Psi)));
        end
        %F(ind) = det(Psi);
        
        
        Psi_(:,:,ind) = Psi;
       
    end

    F = F./Flim;
    %F = F./max(F);
%F = abs(F./max(abs(F)));


end




%% Helper functions

% Leading constant
function c = cij(a,i,j)
    c = exp(-epsilon(a)) / sqrt( prod( 1 + (1:i)/a) * prod( 1 + (1:j)/a) );
end

% Part of leading constant
function e = epsilon(a)
    e = 1/(12*a) - 1/(360*a^3) + 1/(1260*a^6) - 1/(1680*a^7);
end

% exp(phi(a,t))dt is the measure which D_polys are orthogonal w.r.t.
function z = phi(a,t)
    z = a*log1p(t*sqrt(2/a)) - t*sqrt(2*a);
end


% Compute generalized Laguerre coefficients
function L = LaguerreCoef(M,a)
    %A = zeros(M,M);
    L = cell(M,1);
    % Initialize the first two, then do recursion to calculate the rest
    L{1} = 1;
    L{2} = [-1 a+1];
    
    for n=1:M-2
       L{n+2} =  PolyAdd(conv([-1/(n+1), (2*n+1+a) / (n+1)],L{n+1}), -(n+a)/(n+1) * L{n});
    end    
        %A(:,M-k-1) = (2*k+a+1)/(k+1) * A(:,M-k) - (k+a)/(k+1) * A(:,M-k+1) - 1/(k+1)*circshift(A(:,M-k),-1);


end


% This is the "throw it all in the A matrix" version - kind of hard to work
% with because of the indexing. Potentially easier with cells (above
% version)

% function A = LaguerreCoef(M,a)
%     A = zeros(M,M);
%     
%     % Initialize the first two, then do recursion to calculate the rest
%     A(M,M) = 1;
%     A(M,M-1) = a+1;
%     A(M-1,M-1) = -1;
%     
%     for k=1:M-2
%         A(:,M-k-1) = (2*k+a+1)/(k+1) * A(:,M-k) - (k+a)/(k+1) * A(:,M-k+1) - 1/(k+1)*circshift(A(:,M-k),-1);
%     end
% 
% end



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

% Dpolys
%       Computes the coefficients for the D polynomials.
function [D] = Dpolys(a,M)
    q = sqrt(2/a);
	D = cell(1,M);
    D{1} = 1; %D0
    D{2} = [2,-q];
    
    for n=2:M-1
       D{n+1} =  PolyAdd(conv([2, -(2*n+1)*q],D{n}), -(2*(n)+(n)^2*q^2)*D{n-1});
    end
end


% Normalization
%   Computes lim y-> inf F(y) i.e. the normalization constant for the CDF
%   Note that you don't need a very big y for convergence
function [Flim] = Normalization(N, M, mu1, y, D)
    Psi = zeros(M,M);
    
    a = N - M;
    
    %lim = 1.1*max(y);
    lim = 1.03*max(y);
    
    F01 = @( t ) hypergeom( [], N-M+1, mu1*t );
    
        for i = 1:M 
         
            %C_i1 = cij(a,i,1) * a^(-i/2) * a^(-1/2)*sqrt(factorial(i)*factorial(1))/sqrt(pi);
            C_i1 = cij(a,i,1) / sqrt( pi * 2^(i+1) * factorial(i) ); %j=1
            
            for l=i:M      
                H = 0;
                for k=2:M
                    H = H + integral( @(t) exp( phi(a,t) ) .* D{M}(k) .* D{M+1-i}(l-i+1) .* t.^(2*M-l-k), ...
                        -sqrt(a/2), lim );
                end
                
                Psi(i,1) = Psi(i,1) + C_i1 * ( H + integral( @(t) exp( phi(a,t) ) .* D{M}(1) .* D{M+1-i}(l-i+1) .* ...
                    F01( a + t*sqrt(2*a) ) .* t.^(M-l), -sqrt(a/2), lim ) );

            end
                    

            for j=2:M 
               Psi(i,j) = cij(a,i,j) / sqrt( pi * 2^(i+j) * factorial(i) * factorial(j) ) * ...
                    integral( @(t) exp( phi(a,t) ) .* polyval( D{M-i+1}, t ) ...
                    .* polyval( D{M-j+1}, t ), -sqrt(a/2), lim );

            end
        end


    Flim = det(Psi);

end
