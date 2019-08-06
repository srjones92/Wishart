% ## Copyright (C) 2017 Scott Jones <scott.r.jones@asu.edu>
% ##
% ## Derived from wishrnd.m ( (C) 2013 Nir Krakauer <nkrakauer@ccny.cuny.edu>)
% ##
% ## This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
% ##
% ## Octave is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% ##
% ## You should have received a copy of the GNU General Public License along with Octave; see the file COPYING.  If not, see <http://www.gnu.org/licenses/>.
% 
% ## -*- texinfo -*-
% ## @deftypefn  {Function File} {} [@var{W}[, @var{D}]] = wishrnd (@var{Sigma}, @var{N}[, @var{D}][, @var{n_trials}=1])
% ## Return a random matrix sampled from the Wishart distribution with given parameters
% ##
% ## Inputs: the @var{p} x @var{p} positive definite matrix @var{Sigma} and scalar degrees of freedom parameter @var{N} (and optionally the Cholesky factor @var{D} of @var{Sigma}).
% ## @var{N} can be non-integer as long as @var{N} > @var{p}
% ##
% ## Output: a random @var{p} x @var{p}  matrix @var{W} from the Wishart(@var{Sigma}, @var{N}) distribution. If @var{n_trials} > 1, then @var{W} is @var{p} x @var{p} x @var{n_trials} and holds @var{n_trials} such random matrices. (Optionally, the Cholesky factor @var{D} of @var{Sigma} is also returned.)
% ##
% ## Averaged across many samples, the mean of @var{W} should approach @var{N}*@var{Sigma}, and the variance of each element @var{W}_ij should approach @var{N}*(@var{Sigma}_ij^2 + @var{Sigma}_ii*@var{Sigma}_jj)
% ##
% ## Reference: Yu-Cheng Ku and Peter Bloomfield (2010), Generating Random Wishart Matrices with Fractional Degrees of Freedom in OX, http://www.gwu.edu/~forcpgm/YuChengKu-030510final-WishartYu-ChengKu.pdf
% ## 
% ## @seealso{iwishrnd, wishpdf}
% ## @end deftypefn
% 
% ## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
% ## Description: Compute the probability density function of the Wishart distribution

function [W, D] = wishrndC(Sigma, N, D, n_trials)

if (nargin < 2)
  print_usage ();
  n_trials = 1;
end

if nargin < 3 || isempty(D)
  n_trials = 1;
  try
    D = chol(Sigma);
  catch
    error('wishrnd: Cholesky decomposition failed; Sigma probably not positive definite')
  end
end

p = size(D, 1);

if N < p
  N = floor(N); %#distribution not defined for small noninteger N
  N_isint = 1;
else 
%#check for integer degrees of freedom
 N_isint = (N == floor(N));
end

if ~N_isint
  [ii, jj] = ind2sub([p, p], 1:(p*p));
end

if n_trials > 1
  W = nan(p, p, n_trials);
end

for i = 1:n_trials
  if N_isint
	Z = diag(sqrt(gamrnd( N - p+(1:p),1)));

	Z(triu(true(p),1)) = 1/sqrt(2)*(randn(p*(p-1)/2,1)+1i*randn(p*(p-1)/2,1));

	W(:,:,i) = (D*Z)*(D*Z)';

      
      %     Z = 1/sqrt(2)*(randn(N, p) + 1i*randn(N,p)) * D;
%     W(:, :, i) = Z'*Z;
  else
	Z = diag(sqrt(gamrnd( N - p+(1:p),1)));
    
    
	Z(triu(true(p),1)) = 1/sqrt(2)*(randn(p*(p-1)/2,1)+1i*randn(p*(p-1)/2,1));

	W(:,:,i) = (D*Z)*(D*Z)';

%    Z = diag(sqrt(chi2rnd(N - (0:(p-1))))); #fill diagonal
%    #note: chi2rnd(x) is equivalent to 2*randg(x/2), but the latter seems to offer no performance advantage
%    Z(ii > jj) = randn(p*(p-1)/2, 1); #fill lower triangle with normally distributed variates
%    Z = D * Z;
%    W(:, :, i) = Z*Z';
  end

end

end



%!assert(size (wishrnd (1,2,1)), [1, 1]);
%!assert(size (wishrnd ([],2,1)), [1, 1]);
%!assert(size (wishrnd ([3 1; 1 3], 2.00001, [], 1)), [2, 2]);
%!assert(size (wishrnd (eye(2), 2, [], 3)), [2, 2, 3]);

%% Test input validation
%!error wishrnd ()
%!error wishrnd (1)
%!error wishrnd ([1; 1], 2)

