% Log-Normal Distribution for 30% Variance
m = 383.83; % mean
v = 30; % variance
mu = log((m^2)/sqrt(v+m^2))
sigma = sqrt(log(v/(m^2)+1))
% rng('default') % For reproducibility
r = lognrnd(mu,sigma,[1000, 1]);
save -ascii Ae_kTDA.dat r
hist(r)
%%