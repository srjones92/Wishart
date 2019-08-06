%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C_CDF_IS
%
%   Uses importance sampling algorithm from Jiang et al.'s paper
%       "RARE-EVENT ANALYSIS FOR EXTREMAL EIGENVALUES 
%       OF WHITE WISHART MATRICES" 
%   to compute rare event probabilities of the largest eigenvalue of a 
%   central white Wishart matrix. Derived from R code supplied by the 
%   authors of the paper, available at
%   http://users.stat.umn.edu/~xuxxx360/IS.R
%
%   Inputs:      
%       M: size of Wishart matrix / Gaussian covariance parameter
%       N: degrees of freedom
%       x: integral lower limit - we're estimating P(\lambda_1 > N*x)
%       nTrials: specifies number of trials
%
%   Returns:
%       
%
%   Notes:
%       The author's notation has the following correspondence with the  
%       usual seen in this code library:
%           n < P  <-> M < N
%           p = N (degrees of freedom)
%           n = M (size of Wishart matrix M x M)
%           beta: 1 real 2 complex 4 quaternion. Hardcoded to 2 in main
%               function for problems of interest
%   
%       In the initialization, 1 is subtracted from N and M, to construct
%       an M-1 x M-1 "L-matrix" for iterations of the importance sampling.
%
%
%   Original R code included in comments at the end of file, including
%   additional code used for Tracy Widom and direct Monte Carlo estimation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P_l1] = C_CDF_IS(M,N,x, nTrials)

% hard coding it for complex matrices
Beta = 2;
%Beta = 1;

P_l1 = zeros(length(x),1);
std_l1 = zeros(length(x),1);
ratio_l1 = zeros(length(x),1);

[P_l1, std_l1, ratio_l1] = EffSim(Beta, M-1, N-1, x, nTrials);   


% ####   Note here we take n,p as n-1 and p-1
% n <- n-1
% p <- p-1
% iter_num <- 10000 #iteration number of importance sampling estimator
% x <- 1.9 
% 
% EffSim(beta,n,p,x, iter_num) 
% #Importance sampling estimator for P(\lambda_1 > px)
% #Output is (Importance sampling estimator Est., 
% # estimated standard deviations Std., 
% # ratio Std./Est.)
% # Note the standard deviation of the estimate â€œEst.â€ is Std./sqrt(iter_num).





end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   EffSim
%
%   Computes the L estimator given on page 9 of Jiang 2016
%
%   Inputs:
%       Beta: 1 real 2 complex 4 quaternion. Hardcoded to 2 in main
%       function for problems of interest
%       M: size of the L matrix 
%       N: degrees of freedom
%       x:
%       nTrials
%   
%   Returns:
%       LL: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [muLL, stdLL, ratioLL] = EffSim( Beta, M,N, x, nTrials)
% EffSim <- function(beta,n,p,x, iter_num){
% 	LL <- rep(0, iter_num)
% 
% 	for (iter in 1: iter_num){
% 		B <- array(0, c(n,n))
% 
% 		for (i in 1:n){
% 			B[i,i] <- sqrt(rchisq(1,beta*(p+1-i)))
% 		}
% 
% 		for (i in 1:(n-1)){
% 			B[i+1,i] <- sqrt(rchisq(1,beta*(n-i)))
% 		}
% 
% 		L=B%*%t(B)
% 
% 		eigenvalues <- eigen(L)$values
% 		lambda2 <-eigenvalues[1]
% 
% 		lambda1 <- rexp(1,rate=(x-beta)/x/2)+max(lambda2,(p+1)*x)
% 
% 		LL[iter] <- Lfun(1, n+1, p+1, x,lambda1, lambda2, eigenvalues)
% 	}
% 	return(c(mean(LL), 
% 		sqrt(mean(LL^2)-mean(LL)^2),
% 		sqrt(mean(LL^2)-mean(LL)^2)/mean(LL)))
% }
muLL = zeros(size(x));
stdLL = zeros(size(x));

for iter = 1:nTrials
    B = zeros(M,M);
    
    for k = 1:M
        B(k,k) = sqrt(chi2rnd(Beta * (N+1-k)));
        if k < M-1
           B(k+1,k) = sqrt(chi2rnd(Beta * (M-k)));
        end
    end
    
    L = B*B';
    
    eigenvalues = eig(L);
    eigenvalues = sort(eigenvalues,'descend');
    
    lambda2 = eigenvalues(1);
    % MATLAB uses parameter MU which is 1/rate in R
    
    for k=1:length(x)
        lambda1 = exprnd( 2*x(k) / (x(k)-Beta) ) + max(lambda2, (N+1)*x(k) ) ;

        %Lfun(1,M+1,N+1,x,lambda1,lambda2,eigenvalues)
        %LL(iter) = Lfun(Beta,M+1,N+1,x,lambda1,lambda2,eigenvalues);
        Lk(k) = Lfun(Beta,M+1,N+1,x(k),lambda1,lambda2,eigenvalues);
        muLL(k) = muLL(k) + 1/nTrials * Lk(k);
    end
end
     
%     muLL = mean(LL);
%     stdLL = std(LL);
%     ratioLL = stdLL / muLL;
    
ratioLL = stdLL./muLL;

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Lfun
%
%   Computes the L estimator given on page 9 of Jiang 2016
%
%   Inputs:
%       Beta: 1 real 2 complex 4 quaternion. Hardcoded to 2 in main
%       function for problems of interest
%       M: size of the L matrix 
%       N: degrees of freedom
%       lambda1: largest eigenvalue, 
%       lambda2: second largest eigenvalue, computed as eigenvalues(1)
%   
%   Returns:
%       LL: L estimator
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LL] = Lfun( Beta, M, N, x, lambda1, lambda2, eigenvalues)
% Lfun <- function(beta, n, p, x,lambda1,lambda2, eigenvalues){
% 	prob <- 0
% 	prob <- log((x-beta)/x/2)-(x-beta)/x/2*(lambda1-max(lambda2,p*x))-
%log(n)+(beta*(n+p-1)/2)*log(2)-(lgamma(1+beta/2))+(lgamma(1+n*beta/2))+
%(lgamma(p*beta/2))-sum(beta*log(lambda1-eigenvalues))-log(lambda1)*(beta*(p-n+1)/2-1)+lambda1/2
% return(exp(-prob))
% }


% Notes 6/27/18:
% I think this is based on eq 14 on page 9 of the Jiang paper?

% prob = log(M) + log(2) * ( -N*M*Beta/2 + (M-1)*(N-1)*Beta/2 ) + gammaln( 1 + Beta/2 ) ...
%     - gammaln( 1 + M*Beta/2 )  - gammaln( N*Beta/2 ) + Beta * sum( log( lambda1 - eigenvalues ) ) ...
%     + ( Beta*(N-M+1)/2 -1 ) * log(lambda1) - lambda1/2 - log( (x-Beta) / (2*x) ) ...
%     - (x-Beta) / (2*x) * (lambda1 - max(N*x, lambda2));

% ORIGINAL FROM R SCRIPT
prob = log((x-Beta)/x/2) - (x-Beta) / x/2 * (lambda1 - max(lambda2,N*x)) - ...
     log(M) + (Beta*(M+N-1)/2) * log(2) - (gammaln(1+Beta/2)) + (gammaln(1+M*Beta/2)) + ...
     (gammaln(N*Beta/2)) - sum(Beta*log(lambda1-eigenvalues)) - log(lambda1) * ... 
     (Beta*(N-M+1)/2-1) + lambda1/2;

%  if lambda1 < N*x
%      prob = 0;
%  end


 LL = exp(-prob);



end









%%%%%%%%
%
%       ORIGINAL R CODE BELOW
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ###########################################
% # Importance Sampling estimator
% ###########################################
% 
% n <- 10
% p <- 100
% beta <- 1
% 
% Lfun <- function(beta, n, p, x,lambda1,lambda2, eigenvalues){
% 	prob <- 0
% 	prob <- log((x-beta)/x/2)-(x-beta)/x/2*(lambda1-max(lambda2,p*x))-log(n)+(beta*(n+p-1)/2)*log(2)-(lgamma(1+beta/2))+(lgamma(1+n*beta/2))+(lgamma(p*beta/2))-sum(beta*log(lambda1-eigenvalues))-log(lambda1)*(beta*(p-n+1)/2-1)+lambda1/2
% return(exp(-prob))
% }
% 
% 
% EffSim <- function(beta,n,p,x, iter_num){
% 	LL <- rep(0, iter_num)
% 
% 	for (iter in 1: iter_num){
% 		B <- array(0, c(n,n))
% 
% 		for (i in 1:n){
% 			B[i,i] <- sqrt(rchisq(1,beta*(p+1-i)))
% 		}
% 
% 		for (i in 1:(n-1)){
% 			B[i+1,i] <- sqrt(rchisq(1,beta*(n-i)))
% 		}
% 
% 		L=B%*%t(B)
% 
% 		eigenvalues <- eigen(L)$values
% 		lambda2 <-eigenvalues[1]
% 
% 		lambda1 <- rexp(1,rate=(x-beta)/x/2)+max(lambda2,(p+1)*x)
% 
% 		LL[iter] <- Lfun(1, n+1, p+1, x,lambda1, lambda2, eigenvalues)
% 	}
% 	return(c(mean(LL), 
% 		sqrt(mean(LL^2)-mean(LL)^2),
% 		sqrt(mean(LL^2)-mean(LL)^2)/mean(LL)))
% }
% 
% ####   Note here we take n,p as n-1 and p-1
% n <- n-1
% p <- p-1
% iter_num <- 10000 #iteration number of importance sampling estimator
% x <- 1.9 
% 
% EffSim(beta,n,p,x, iter_num) 
% #Importance sampling estimator for P(\lambda_1 > px)
% #Output is (Importance sampling estimator Est., 
% # estimated standard deviations Std., 
% # ratio Std./Est.)
% # Note the standard deviation of the estimate â€œEst.â€ is Std./sqrt(iter_num).
% 
% 
% ###########################################
% # T-W estimator
% ###########################################
% 
% library(RMTstat)
% 
% n <- 10
% p <- 100
% beta <- 1
% x<-1.9
% 
% # TW approximation has been proposed in Johnstone (2001)
% bb <- (p*x-(sqrt(n)+sqrt(p-1))^2)/(sqrt(n)+sqrt(p-1))/(1/sqrt(n)+1/sqrt(p-1))^{1/3}
% 1-ptw(bb)
% 
% # TW approximation has been proposed in Ma (2012), 
% bb <- (p*x-(sqrt(n-.5)+sqrt(p-.5))^2)/(sqrt(n-.5)+sqrt(p-.5))/(1/sqrt(n-.5)+1/sqrt(p-.5))^{1/3}
% 1-ptw(bb)
% 
% 
% 
% ###########################################
% # Direct MC estimator 
% ###########################################
% 
% n <- 10
% p <- 100
% beta <- 1
% 
% iter_num <- 1000000
% lambda1 <- rep(0, iter_num)
% 
% for (iter in 1: iter_num){
% B <- array(0, c(n,n))
% 
% for (i in 1:n){
% 	B[i,i] <- sqrt(rchisq(1,beta*(p+1-i)))
% }
% 
% for (i in 1:(n-1)){
% 	B[i+1,i] <- sqrt(rchisq(1,beta*(n-i)))
% }
% 
% L=B%*%t(B)
% 
% eigenvalues <- eigen(L)$values
% lambda1[iter] <-eigenvalues[1]
% }
% 
% x <- 1.9
% sum(lambda1/p > x)/iter_num
