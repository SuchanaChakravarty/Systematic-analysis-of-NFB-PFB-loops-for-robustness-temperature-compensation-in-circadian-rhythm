function dcdt = osc_Luca(t,c,Temp,...
    Ae_p0,Ae_p3,Ae_d1,...
    Ae_d2,Ae_k01,Ae_k02,Ae_k03,...
    Ae_k04)
% variables
OO = c(1);
PP = c(2);
OP = c(3);
PO = c(4);
B = c(5);
A = c(6);
Tot_Rate_Var = c(7);
% parameters, random pre-exponential factor(Ae)
% Generate Random Number first before Running this code
p0 = Ae_p0*exp(-(6.45*10^3)/(8.3144598*Temp)); % 28.42
p3 = Ae_p3*exp(-(10*10^3)/(8.3144598*Temp)); %6.782 
d1 = Ae_d1*exp(-(10*10^3)/(8.3144598*Temp)); %6.782
d2 = Ae_d2*exp(-(6.45*10^3)/(8.3144598*Temp)); %28.42 
k01 = Ae_k01*exp(-(21.55*10^3)/(8.3144598*Temp));  %0.06714
k02 = Ae_k02*exp(-(21.55*10^3)/(8.3144598*Temp));  %0.06714
k03 = Ae_k03*exp(-(21.55*10^3)/(8.3144598*Temp));  %0.06714
k04 = Ae_k04*exp(-(21.55*10^3)/(8.3144598*Temp));  %0.06714
% parameter as base value at each temperature, pre-exponential factor
% =Ae=383.83
p0p = 383.83*exp(-(6.45*10^3)/(8.3144598*Temp)); % 28.42
p3p = 383.83*exp(-(10*10^3)/(8.3144598*Temp)); %6.782 
d1p = 383.83*exp(-(10*10^3)/(8.3144598*Temp)); %6.782
d2p = 383.83*exp(-(6.45*10^3)/(8.3144598*Temp)); %28.42 
k01p = 383.83*exp(-(21.55*10^3)/(8.3144598*Temp));  %0.06714
k02p = 383.83*exp(-(21.55*10^3)/(8.3144598*Temp));  %0.06714
k03p = 383.83*exp(-(21.55*10^3)/(8.3144598*Temp));  %0.06714
k04p = 383.83*exp(-(21.55*10^3)/(8.3144598*Temp));  %0.06714
%total no of parameter
n=8;
%total parameter variation
Z=((abs(p0-p0p)/p0p)+(abs(p3-p3p)/p3p)+(abs(d1-d1p)/d1p)...
    +(abs(d2-d2p)/d2p)+(abs(k01-k01p)/k01p)+(abs(k02-k02p)/k02p)...
    +(abs(k03-k03p)/k03p)+(abs(k04-k04p)/k04p))/n;
% initialization 
dcdt = zeros(7,1);
% differential equations
dcdt(1) = (k02.*B.*PP + d1.*OO.*PO - k01.*A.*OO - p0.*OO.*PP);
dcdt(2) = (k01.*A.*OO + p3.*OP.*PP - k02.*B.*PP - d2.*OO.*PP);
dcdt(3) = (p0.*OO.*PP - p3.*OP.*PP);
dcdt(4) = (d2.*OO.*PP - d1.*OO.*PO);
dcdt(5) = (k03.*A.*PP - k04.*B.*OO);
dcdt(6) = (k04.*B.*OO - k03.*A.*PP);
dcdt(7) = (Z- Tot_Rate_Var);