t0 = 0;
tf = 1150;
tspan = t0:0.02:tf;
css = [0.48 0.0075]';
Temp = zeros(7,1);
Y_val = zeros(length(tspan),7);
U_val = zeros(length(tspan),7);
for i=1:7
    i
    Temp(i,1)=278+5*i;
[t,c]=ode45(@(t,c) subdep(t,c,Temp(i,1)) ,tspan,css);
% figure(2)
% plot(t,c(:,1),'.-','Color',[rand,rand,rand],'DisplayName','Y')
% hold on
% plot(t,c(:,2),'.-','Color',[rand,rand,rand],'DisplayName','U')
Y_val(:,i) = c(:,1);
U_val(:,i) = c(:,2);
[peakval,locval]=findpeaks(U_val(:,i),t);
period(:,i) = mean(diff(locval));
period1 = period';
end
data = [Temp period1];
save -ascii period_subdep_k1_k2_fixed.dat data
%%
