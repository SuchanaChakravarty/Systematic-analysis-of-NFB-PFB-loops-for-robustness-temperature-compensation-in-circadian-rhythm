function dcdt = osc_Luca(t,c,Temp)
% variables
OO = c(1);
PP = c(2);
OP = c(3);
PO = c(4);
B = c(5);
A = c(6);
Temp1 = 298; % the rates that are temp compensated kept fixed at 298k
% parameters
p0 = 383.83*exp(-(6.45*10^3)/(8.3144598*Temp)); % 28.42
p3 = 383.83*exp(-(10*10^3)/(8.3144598*Temp)); %6.782 
d1 = 383.83*exp(-(10*10^3)/(8.3144598*Temp1)); %6.782
d2 = 383.83*exp(-(6.45*10^3)/(8.3144598*Temp1)); %28.42 
k01 = 383.83*exp(-(21.55*10^3)/(8.3144598*Temp));  %0.06714
k02 = 383.83*exp(-(21.55*10^3)/(8.3144598*Temp));  %0.06714
k03 = 383.83*exp(-(21.55*10^3)/(8.3144598*Temp));  %0.06714
k04 = 383.83*exp(-(21.55*10^3)/(8.3144598*Temp));  %0.06714
% initialization 
dcdt = zeros(6,1);
% differential equations
dcdt(1) = (k02.*B.*PP + d1.*OO.*PO - k01.*A.*OO - p0.*OO.*PP);
dcdt(2) = (k01.*A.*OO + p3.*OP.*PP - k02.*B.*PP - d2.*OO.*PP);
dcdt(3) = (p0.*OO.*PP - p3.*OP.*PP);
dcdt(4) = (d2.*OO.*PP - d1.*OO.*PO);
dcdt(5) = (k03.*A.*PP - k04.*B.*OO);
dcdt(6) = (k04.*B.*OO - k03.*A.*PP);