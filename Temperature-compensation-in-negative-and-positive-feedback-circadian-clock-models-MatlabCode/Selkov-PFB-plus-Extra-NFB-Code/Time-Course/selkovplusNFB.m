function dcdt = selkovplusNFB(t,c,Temp)
% variables
Y = c(1);
U = c(2);
% parameters : Change K3, K5 and K1 values according to the cases
% k1 = 383.83.*exp(-(20.877.*10^3)/(8.3144598.*Temp)); % 0.0841 (k3 k5 both Large)
% k1 = 383.83.*exp(-(16.35.*10^3)/(8.3144598.*Temp)); % 0.2325 (Large k5 Small k3)
k1 = 383.83.*exp(-(21.5.*10^3)/(8.3144598.*Temp)); % 0.0691 (Large k3 Small k5)
k2 = 383.83.*exp(-(27.5.*10^3)/(8.3144598.*Temp)); % 0.0058 
k3 = 383.83.*exp(-(5.17*10^3)/(8.3144598.*Temp)); % 47.6363 (Large K3)
% k3 =383.83.*exp(-(12.0203*10^3)/(8.3144598.*Temp)); % 3 (Small K3);
k4 = 383.83.*exp(-(14.8691.*10^3)/(8.3144598.*Temp)); % 0.9503 
% k5 = 383.83.*exp(-(5.17*10^3)/(8.3144598.*Temp)); % 47.6363 (Large K5)
k5 =383.83.*exp(-(12.0203*10^3)/(8.3144598.*Temp)); % 3 (Small K5);

% initialization 
dcdt = zeros(2,1);
% differential equations
dcdt(1) = ((k1./(1+(k5.*U))) - k2.*Y - k3.*U^2.*Y);
dcdt(2) = (k2.*Y + k3.*U^2.*Y - k4.*U);