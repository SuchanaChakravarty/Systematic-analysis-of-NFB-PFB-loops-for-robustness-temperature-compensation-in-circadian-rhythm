t0 = 0;
tf = 1150;
tspan = t0:0.02:tf;
css = [0.68 1.36 0.34 0.0]'; 
load 'Ae_khalf.dat';
load 'Ae_kUTA.dat';
load 'Ae_kDTA.dat';
load 'Ae_kTU0.dat';
load 'Ae_kTUA.dat';
load 'Ae_kTDA.dat';
load 'Ae_kSDA.dat';
load 'Ae_kDS0.dat';
load 'Ae_kDSA.dat';
load 'Ae_kUSA.dat';
load 'Ae_kSU0.dat';
load 'Ae_kSUA.dat';
Tot_Rate_Var_val = zeros(7,length(Ae_khalf));
for q=1:length(Ae_khalf)
    q
Ae_khalf(q,1)=Ae_khalf(q);
Ae_kUTA(q,1)=Ae_kUTA(q);
Ae_kDTA(q,1)=Ae_kDTA(q);
Ae_kTU0(q,1)=Ae_kTU0(q);
Ae_kTUA(q,1)=Ae_kTUA(q);
Ae_kTDA(q,1)=Ae_kTDA(q);
Ae_kSDA(q,1)=Ae_kSDA(q);
Ae_kDS0(q,1)=Ae_kDS0(q);
Ae_kDSA(q,1)=Ae_kDSA(q);
Ae_kUSA(q,1)=Ae_kUSA(q);
Ae_kSU0(q,1)=Ae_kSU0(q);
Ae_kSUA(q,1)=Ae_kSUA(q);
Temp = zeros(7,1);
T_val = zeros(length(tspan),7);
ST_val = zeros(length(tspan),7);
S_val = zeros(length(tspan),7);
M_val = zeros(length(tspan),7);
for i=1:7
    i
    Temp(i,1)=278+5*i;
[t,c]=ode45(@(t,c) kaiABC_Rust(t,c,Temp(i,1),Ae_khalf(q,1),...
    Ae_kUTA(q,1),Ae_kDTA(q,1),...
    Ae_kTU0(q,1),Ae_kTUA(q,1),Ae_kTDA(q,1),Ae_kSDA(q,1),...
    Ae_kDS0(q,1),Ae_kDSA(q,1),Ae_kUSA(q,1),Ae_kSU0(q,1),...
    Ae_kSUA(q,1)),tspan,css);
% figure(2)
% plot(t,c(:,1),'.-','Color',[rand,rand,rand],'DisplayName','T')
% hold on
% plot(t,c(:,2),'.-','Color',[rand,rand,rand],'DisplayName','ST')
% plot(t,c(:,3),'.-','Color',[rand,rand,rand],'DisplayName','S')
T_val(:,i) = c(:,1);
ST_val(:,i) = c(:,2);
S_val(:,i) = c(:,3);
M_val(:,i) = c(28750,4);
[peakval,locval]=findpeaks(ST_val(:,i),t);
period(:,i) = mean(diff(locval));
end
period1(:,q) = period;
M_val1(:,q) = M_val((length(tspan)-1)/2,1:7);
Tot_Rate_Var_val(:,q) =M_val1(:,q);
end
save -ascii Rust_period.dat period1
save -ascii Tot_Rate_Var_rust.dat Tot_Rate_Var_val
%%
