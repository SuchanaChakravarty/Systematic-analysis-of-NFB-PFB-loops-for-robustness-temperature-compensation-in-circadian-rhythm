function dcdt = NFL(t,c,Temp,...
    Ae_k1,Ae_k2,Ae_d1,...
    Ae_d2,Ae_k,Ae_a1,Ae_a2)
% variables
X = c(1);
Y = c(2);
Tot_Rate_Var = c(3);
% parameters, random pre-exponential factor(Ae)
% Generate Random Number first before Running this code
k1 = Ae_k1.*exp(-(31.8565.*10^3)/(8.3144598.*Temp)); %k1=0.001
k2 = Ae_k2.*exp(-(31.8565.*10^3)/(8.3144598.*Temp)); %k2=0.001
d1 = Ae_d1.*exp(-(14.7420.*10^3)/(8.3144598.*Temp)); %d1=1.0003
d2 = Ae_d2.*exp(-(14.7420.*10^3)/(8.3144598.*Temp)); %d2=1.0003
k = Ae_k.*exp(-(12.*10^3)/(8.3144598.*Temp)); %k=3.0253
a1 = Ae_a1.*exp(-(14.01.*10^3)/(8.3144598.*Temp)); %a1=1.3442
a2 = Ae_a2.*exp(-(17.9420.*10^3)/(8.3144598.*Temp)); %a2=0.2749
% parameter as base value at each temperature, pre-exponential factor
% =Ae=383.83
k1_p = 383.83.*exp(-(31.8565.*10^3)/(8.3144598.*Temp)); 
k2_p = 383.83.*exp(-(31.8565.*10^3)/(8.3144598.*Temp));
d1_p = 383.83.*exp(-(14.7420.*10^3)/(8.3144598.*Temp)); 
d2_p = 383.83.*exp(-(14.7420.*10^3)/(8.3144598.*Temp)); 
k_p = 383.83.*exp(-(12.*10^3)/(8.3144598.*Temp)); 
a1_p = 383.83.*exp(-(14.01.*10^3)/(8.3144598.*Temp)); 
a2_p = 383.83.*exp(-(17.9420.*10^3)/(8.3144598.*Temp));
% total no of parameter
n=7;
%total parameter variation
Z=((abs(k1-k1_p)/k1_p)+(abs(k2-k2_p)/k2_p)+(abs(d1-d1_p)/d1_p)...
    +(abs(d2-d2_p)/d2_p)+(abs(k-k_p)/k_p)+(abs(a1-a1_p)/a1_p)...
    +(abs(a2-a2_p)/a2_p))/n;
% initialization 
dcdt = zeros(3,1);
% differential equations
dcdt(1) = (((a1.*k)./(k+Y)) - ((d1.*X)./(k1+X)));
dcdt(2) = ((a2.*X) - ((d2.*Y)./(k2+Y)));
dcdt(3) = (Z- Tot_Rate_Var);