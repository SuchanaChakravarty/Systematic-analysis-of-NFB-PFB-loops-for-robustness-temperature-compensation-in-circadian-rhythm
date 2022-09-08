t0 = 0;
tf = 1150;
tspan = t0:0.02:tf;
css = [0.68 1.36 0.34]';
Temp = zeros(7,1);
T_val = zeros(length(tspan),7);
ST_val = zeros(length(tspan),7);
S_val = zeros(length(tspan),7);
for i=1:7
    i
    Temp(i,1)=278+5*i;
[t,c]=ode45(@(t,c) kaiABC_Rust(t,c,Temp(i,1)),tspan,css);
% figure(2)
% plot(t,c(:,1),'.-','Color',[rand,rand,rand],'DisplayName','T')
% hold on
% plot(t,c(:,2),'.-','Color',[rand,rand,rand],'DisplayName','ST')
% plot(t,c(:,3),'.-','Color',[rand,rand,rand],'DisplayName','S')
T_val(:,i) = c(:,1);
ST_val(:,i) = c(:,2);
S_val(:,i) = c(:,3);
[peakval,locval]=findpeaks(ST_val(:,i),t);
period(:,i) = max(diff(locval));
period1 = period';
end
data = [Temp period1];
save -ascii period_kaiABC_Rust_kUTA_fixed.dat data
%%