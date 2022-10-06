t0 = 0;
tf = 1150;
tspan = t0:0.02:tf;
Temp = zeros(7,1);
OO_val = zeros(length(tspan),7);
PP_val = zeros(length(tspan),7);
OP_val = zeros(length(tspan),7);
PO_val = zeros(length(tspan),7);
B_val = zeros(length(tspan),7);
A_val = zeros(length(tspan),7);
css = [2 1 0 0 2 1]';
for i=1:7
    i
    Temp(i,1)=278+5*i;
[t,c]=ode45(@(t,c) osc_Luca(t,c,Temp(i,1)) ,tspan,css);
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
[peakval,locval]=findpeaks(PP_val(:,i),t);
period(:,i) = max(diff(locval));
period1 = period';
end
data = [Temp period1];
save -ascii period_osc_Luca_d1_d2_fixed.dat data
%%