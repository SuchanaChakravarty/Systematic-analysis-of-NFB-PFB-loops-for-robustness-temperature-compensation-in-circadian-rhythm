function dcdt = NFL(t,c,Temp)
% variables
X = c(1);
Y = c(2);
% parameters
k1 = 383.83.*exp(-(31.8565.*10^3)/(8.3144598.*Temp)); %k1=0.001
k2 = 383.83.*exp(-(31.8565.*10^3)/(8.3144598.*Temp)); %k2=0.001
d1 = 383.83.*exp(-(14.7420.*10^3)/(8.3144598.*Temp)); %d1=1.0003
d2 = 383.83.*exp(-(14.7420.*10^3)/(8.3144598.*Temp)); %d2=1.0003
k = 383.83.*exp(-(12.*10^3)/(8.3144598.*Temp)); %k=3.0253
a1 = 383.83.*exp(-(14.01.*10^3)/(8.3144598.*Temp)); %a1=1.3442
a2 = 383.83.*exp(-(17.9420.*10^3)/(8.3144598.*Temp)); %a2=0.2749
% initialization 
dcdt = zeros(2,1);
% differential equations
dcdt(1) = (((a1.*k)./(k+Y)) - ((d1.*X)./(k1+X)));
dcdt(2) = ((a2.*X) - ((d2.*Y)./(k2+Y)));