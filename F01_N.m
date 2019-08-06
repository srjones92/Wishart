
function z = F01_N(mu1,a,t)

% want to compute 0F1(a+1,mu1(a+t*sqrt(2*a)) for large a and t - this
% introduces a normalization term that is cancelled in the
% numerator/denominator when taking determinants to compute lambda1

max_iter = 1000;

z = 0;

order = 1.3*(ceil(log10(a)) + ceil(log10(mu1)))^2;
for k=0:max_iter

    
    if log(a+t*sqrt(2*a)) > -inf
        A = k*(log(mu1) + log(a+t*sqrt(2*a))) + gammaln(a+1) - gammaln(a+1+k) - gammaln(k+1) - order * log(mu1) ;
    else
        A = - order * log(mu1) ; 
    end
    
    z = z + exp(A);
    
%     if (pochhammer(a+1,k)*factorial(k) < inf) && ( mu1^(k-order^2) * (a+t*sqrt(2*a))^(k) < inf )
%         z = z + mu1^(k-order^2) * (a+t*sqrt(2*a))^(k) / ( pochhammer(a+1,k) * factorial(k));   
%     else
%         break;
%     end
%     
    %z = z + mu1^(k) * (a+t*sqrt(2*a))^(k) / ( prod( a+1+[1:k]) * factorial(k));
end



end