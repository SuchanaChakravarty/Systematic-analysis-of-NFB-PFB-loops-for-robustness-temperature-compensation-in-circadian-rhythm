t0 = 0;
tf = 1150;
tspan = t0:0.02:tf;
css = [0.48 0.0075 0]';
load 'Ae_k1.dat';
load 'Ae_k2.dat';
load 'Ae_k3.dat';
load 'Ae_k4.dat';
load 'Ae_k5.dat';
Tot_Rate_Var_val = zeros(7,length(Ae_k1));
for q=1:length(Ae_k1)
    q
Ae_k1(q,1)=Ae_k1(q);
Ae_k2(q,1)=Ae_k2(q);
Ae_k3(q,1)=Ae_k3(q);
Ae_k4(q,1)=Ae_k4(q);
Ae_k5(q,1)=Ae_k5(q);
Temp = zeros(7,1);
Y_val = zeros(length(tspan),7);
U_val = zeros(length(tspan),7);
M_val = zeros(length(tspan),7);
for i=1:7
    i
    Temp(i,1)=278+5*i;
[t,c]=ode45(@(t,c) selkovplusNFB(t,c,Temp(i,1),Ae_k1(q,1),...
    Ae_k2(q,1),Ae_k3(q,1),Ae_k4(q,1),Ae_k5(q,1)) ,tspan,css);
% figure(2)
% plot(t,c(:,1),'.-','Color',[rand,rand,rand],'DisplayName','Y')
% hold on
plot(t,c(:,2),'.-','Color',[rand,rand,rand],'DisplayName','U')
Y_val(:,i) = c(:,1);
U_val(:,i) = c(:,2);
M_val(:,i) = c(28750,3);
[peakval,locval]=findpeaks(U_val(:,i),t);
period(:,i) = mean(diff(locval));
end
period1(:,q) = period;
M_val1(:,q) = M_val((length(tspan)-1)/2,1:7);
Tot_Rate_Var_val(:,q) =M_val1(:,q);
end
% save -ascii subdep_period_lk3_sk5.dat period1
% save -ascii Tot_Rate_Var_selkov_lk3_sk5.dat Tot_Rate_Var_val
%%
