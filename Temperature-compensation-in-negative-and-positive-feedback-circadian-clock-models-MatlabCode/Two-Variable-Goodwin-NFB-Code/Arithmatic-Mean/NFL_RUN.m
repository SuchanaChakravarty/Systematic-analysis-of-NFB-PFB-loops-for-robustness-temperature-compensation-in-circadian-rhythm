t0 = 0;
tf = 1150;
tspan = t0:0.02:tf;
css = [1 0 0]';
load 'Ae_k1.dat';
load 'Ae_k2.dat';
load 'Ae_d1.dat';
load 'Ae_d2.dat';
load 'Ae_k.dat';
load 'Ae_a1.dat';
load 'Ae_a2.dat';
Tot_Rate_Var_val = zeros(7,length(Ae_k));
for q=1:length(Ae_k)
    q
Ae_k1(q,1)=Ae_k1(q);
Ae_k2(q,1)=Ae_k2(q);
Ae_d1(q,1)=Ae_d1(q);
Ae_d2(q,1)=Ae_d2(q);
Ae_k(q,1)=Ae_k(q);
Ae_a1(q,1)=Ae_a1(q);
Ae_a2(q,1)=Ae_a2(q);
Temp = zeros(7,1);
X_val = zeros(length(tspan),7);
Y_val = zeros(length(tspan),7);
M_val = zeros(length(tspan),7);
for i=1:7
    i
    Temp(i,1)=278+5*i;
[t,c]=ode45(@(t,c) NFL(t,c,Temp(i,1),Ae_k1(q,1),...
    Ae_k2(q,1),Ae_d1(q,1),Ae_d2(q,1),...
    Ae_k(q,1),Ae_a1(q,1),Ae_a2(q,1)),tspan,css);
% figure(2)
% plot(t,c(:,1),'.-','Color',[rand,rand,rand],'DisplayName','X')
% hold on
% plot(t,c(:,2),'.-','Color',[rand,rand,rand],'DisplayName','Y')
X_val(:,i) = c(:,1);
Y_val(:,i) = c(:,2);
M_val(:,i) = c(20000,3);
[peakval,locval]=findpeaks(X_val(:,i),t);
period(:,i) = max(diff(locval));
end
period1(:,q) = period;
M_val1(:,q) = M_val((length(tspan)-1)/2,1:7);
Tot_Rate_Var_val(:,q) =M_val1(:,q);
end
save -ascii small_NFL_period.dat period1
save -ascii Tot_Rate_Var_small_NFL.dat Tot_Rate_Var_val
%%