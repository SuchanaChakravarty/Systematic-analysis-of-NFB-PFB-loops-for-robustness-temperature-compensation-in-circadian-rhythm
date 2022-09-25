function dcdt = selkovplusNFB(t,c,Temp,...
    Ae_k1,Ae_k2,Ae_k3,...
    Ae_k4,Ae_k5)
% variables
Y = c(1);
U = c(2);
Tot_Rate_Var = c(3);
% parameters, random pre-exponential factor(Ae)
% Generate Random Number first before Running this code
k1 = Ae_k1.*exp(-(16.35.*10^3)/(8.3144598.*Temp)); % 0.2325 (Large k5 Small k3)
k2 = Ae_k2.*exp(-(27.5.*10^3)/(8.3144598.*Temp)); % 0.0058 
k3 = Ae_k3.*exp(-(12.0203*10^3)/(8.3144598.*Temp)); % 3 (Small k3);
k4 = Ae_k4.*exp(-(14.8691.*10^3)/(8.3144598.*Temp)); % 0.9503 
k5 = Ae_k5.*exp(-(5.17*10^3)/(8.3144598.*Temp)); % 47.6363 (Large k5)

% parameter as base value at each temperature, pre-exponential factor
% =Ae=383.83
k1p = 383.83.*exp(-(16.35.*10^3)/(8.3144598.*Temp)); % 0.2325 (Large k5 Small k3)
k2p = 383.83.*exp(-(27.5.*10^3)/(8.3144598.*Temp)); % 0.0058 
k3p =383.83.*exp(-(12.0203*10^3)/(8.3144598.*Temp)); % 3 (Small k3);
k4p = 383.83.*exp(-(14.8691.*10^3)/(8.3144598.*Temp)); % 0.9503 
k5p = 383.83.*exp(-(5.17*10^3)/(8.3144598.*Temp)); % 47.6363 (Large k5)
%total no of parameter
n=5;
%total parameter variation
Z=((abs(k1-k1p)/k1p)+(abs(k2-k2p)/k2p)+(abs(k3-k3p)/k3p)...
    +(abs(k4-k4p)/k4p)+(abs(k5-k5p)/k5p))/n;
% initialization 
dcdt = zeros(3,1);
% differential equations
dcdt(1) = ((k1./(1+(k5.*U))) - k2.*Y - k3.*U^2.*Y);
dcdt(2) = (k2.*Y + k3.*U^2.*Y - k4.*U);
dcdt(3) = (Z- Tot_Rate_Var);