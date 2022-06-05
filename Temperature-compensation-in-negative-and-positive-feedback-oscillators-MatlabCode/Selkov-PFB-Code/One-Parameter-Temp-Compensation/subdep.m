function dcdt = subdep(t,c,Temp)
% variables
Y = c(1);
U = c(2);
Temp1 = 298; % the rate that is temp compensated kept fixed at 298k
% parameters
k1 = 383.83.*exp(-(21.5.*10^3)/(8.3144598.*Temp1)); % 0.0654
k2 = 383.83.*exp(-(27.5.*10^3)/(8.3144598.*Temp)); % 0.0058
k3 = 383.83.*exp(-(5.1768.*10^3)/(8.3144598.*Temp)); % 47.5057
k4 = 383.83.*exp(-(14.8691.*10^3)/(8.3144598.*Temp)); % 0.9503
% initialization 
dcdt = zeros(2,1);
% differential equations
dcdt(1) = (k1 - k2.*Y - k3.*U^2.*Y);
dcdt(2) = (k2.*Y + k3.*U^2.*Y - k4.*U);