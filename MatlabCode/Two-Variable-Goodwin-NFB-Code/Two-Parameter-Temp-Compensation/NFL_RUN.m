t0 = 0;
tf = 1150;
tspan = t0:0.02:tf;
css = [1 0]';
Temp = zeros(7,1);
X_val = zeros(length(tspan),7);
Y_val = zeros(length(tspan),7);
for i=1:7
    i
    Temp(i,1)=278+5*i;
[t,c]=ode45(@(t,c) NFL(t,c,Temp(i,1)),tspan,css);
figure(2)
plot(t,c(:,1),'.-','Color',[rand,rand,rand],'DisplayName','X')
hold on
plot(t,c(:,2),'.-','Color',[rand,rand,rand],'DisplayName','Y')
X_val(:,i) = c(:,1);
Y_val(:,i) = c(:,2);
[peakval,locval]=findpeaks(X_val(:,i),t);
period(:,i) = mean(diff(locval));
period1 = period';
end
data = [Temp period1];
save -ascii NFL_period_a1_a2_fixed.dat data
%%
