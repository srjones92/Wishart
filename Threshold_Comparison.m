%% Threshold Comparison
% Compares Monte Carlo, D polynomial analytical, and Tracy Widom
% approximation thresholds


%% Parameters

N = 1e5;
M = 4;

Pf = logspace(-6,-4,3);


%% Monte Carlo Simulation

nTrials = 1e8;
lambda1 = zeros(nTrials,1);
for k=1:nTrials
    lambda1(k) = max(eig(wishrndC(eye(M),N,eye(M),1)));    
end

[FMC, XMC] = ecdf(lambda1);
PFMC = 1 - FMC;

for l=1:length(Pf)
[~, ind(l)] = min(abs(PFMC - Pf(l)));
end

T_MC = XMC(ind);



%% D Polynomial Thresholds

T_D = D_Threshold(M,N,Pf);



%% Tracy Widom Approximation
mu = (sqrt(N)+sqrt(M))^2;
sigma = sqrt(mu)*(1/sqrt(N) + 1/sqrt(M))^(1/3);

Tx = linspace(0.9*min(T_MC), 1.1*max(T_MC), 200);

% shift into corect range. ptw also requires values above and below certain
% indices in the lookup table for some reason; results from some design
% decision in the original R code
Ty = [-17.7971, (Tx - mu)/sigma, 40]; 

% beta = 2 -> Complex
p_TW = ptw(Ty,2);
Pf_TW = 1 - p_TW(2:end-1); % can throw away the end values now

for k=1:length(Pf)
[~, ind_TW(k)] = min(abs(Pf(k)-Pf_TW));
end
T_TW = Tx(ind_TW)';


%% Generate Plots

figure(1);

plot(T_MC,T_MC,'-x',T_MC,T_D,'-*',T_MC,T_TW,'-+');
legend('Monte Carlo','D Polynomial','Tracy-Widom')
%title('Comparison of Threshold Computation Methods with Monte Carlo Experimental Thresholds');
xlabel('Monte Carlo Threshold')
ylabel('Computed Threshold');


figure(2); 

loglog(Pf,abs(T_MC - T_D)./T_MC,'-*',Pf,abs(T_MC-T_TW)./T_MC,'-+')
xlabel('P_F')
ylabel('Threshold Relative Error')
legend('D Polynomial','Tracy-Widom')


%% Check P_F given by plugging in T_MC into other methods

PFMC = 1-FMC(ind);
PFD = 1-C_CDF_D(M,N,(T_MC-a)/sqrt(2*a));
PFTW = 1 - ptw( [-17.791, (T_TW' - mu)/sigma, 40],2);
PFTW = PFTW(2:end-1);


figure(3);
loglog(Pf,PFMC,'-x',Pf,PFD,'-*',Pf,PFTW,'-+');
ylim([1e-6,1e-1])
legend('Monte Carlo', 'D Polynomial','Tracy-Widom');
xlabel('P_F');
ylabel('Computed P_F');


figure(4);
loglog(Pf,abs(Pf-PFD')./Pf, '-*', Pf,abs(Pf-PFTW')./Pf,'-+');
legend('D Polynomial','Tracy-Widom');
xlabel('P_F');
ylabel('P_F Relative Error');
