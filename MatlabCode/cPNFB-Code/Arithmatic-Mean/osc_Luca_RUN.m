t0 = 0;
tf = 1150;
tspan = t0:0.02:tf;
css = [2 1 0 0 2 1 0]';
load 'Ae_p0.dat';
load 'Ae_p3.dat';
load 'Ae_d1.dat';
load 'Ae_d2.dat';
load 'Ae_k01.dat';
load 'Ae_k02.dat';
load 'Ae_k03.dat';
load 'Ae_k04.dat';
Tot_Rate_Var_val = zeros(7,length(Ae_k01));
for q=1:length(Ae_k01)
    q
Ae_p0(q,1)=Ae_p0(q);
Ae_p3(q,1)=Ae_p3(q);
Ae_d1(q,1)=Ae_d1(q);
Ae_d2(q,1)=Ae_d2(q);
Ae_k01(q,1)=Ae_k01(q);
Ae_k02(q,1)=Ae_k02(q);
Ae_k03(q,1)=Ae_k03(q);
Ae_k04(q,1)=Ae_k04(q);
Temp = zeros(7,1);
OO_val = zeros(length(tspan),7);
PP_val = zeros(length(tspan),7);
OP_val = zeros(length(tspan),7);
PO_val = zeros(length(tspan),7);
B_val = zeros(length(tspan),7);
A_val = zeros(length(tspan),7);
M_val = zeros(length(tspan),7);
for i=1:7
    i
    Temp(i,1)=278+5*i;
[t,c]=ode45(@(t,c) osc_Luca(t,c,Temp(i,1),Ae_p0(q,1),...
    Ae_p3(q,1),Ae_d1(q,1),Ae_d2(q,1),...
    Ae_k01(q,1),Ae_k02(q,1),Ae_k03(q,1),...
    Ae_k04(q,1)) ,tspan,css);
% figure(2)
% plot(t,c(:,1),'.-','Color',[rand,rand,rand],'DisplayName','OO')
% hold on
% plot(t,c(:,2),'.-','Color',[rand,rand,rand],'DisplayName','PP')
% plot(t,c(:,3),'.-','Color',[rand,rand,rand],'DisplayName','OP')
% plot(t,c(:,4),'.-','Color',[rand,rand,rand],'DisplayName','PO')
% plot(t,c(:,5),'.-','Color',[rand,rand,rand],'DisplayName','B')
% plot(t,c(:,6),'.-','Color',[rand,rand,rand],'DisplayName','A')
OO_val(:,i) = c(:,1);
PP_val(:,i) = c(:,2);
OP_val(:,i) = c(:,3);
PO_val(:,i) = c(:,4);
B_val(:,i) = c(:,5);
A_val(:,i) = c(:,6);
M_val(:,i) = c(28750,7);
[peakval,locval]=findpeaks(PP_val(:,i),t);
period(:,i) = mean(diff(locval));
end
period1(:,q) = period;
M_val1(:,q) = M_val((length(tspan)-1)/2,1:7);
Tot_Rate_Var_val(:,q) =M_val1(:,q);
end
save -ascii Luca_osc_period.dat period1
save -ascii Tot_Rate_Var_luca.dat Tot_Rate_Var_val
%%
