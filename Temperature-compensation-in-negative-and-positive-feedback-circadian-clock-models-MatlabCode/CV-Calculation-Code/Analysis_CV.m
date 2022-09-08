load Luca_osc_period.dat; % Obtain after Arithmetic Mean Calculation
load small_NFL_period.dat; % Obtain after Arithmetic Mean Calculation
load Rust_period.dat; % Obtain after Arithmetic Mean Calculation
load subdep_period.dat; % Obtain after Arithmetic Mean Calculation
load Tot_Rate_Var_luca.dat; % Obtain after Arithmetic Mean Calculation
load Tot_Rate_Var_small_NFL.dat; % Obtain after Arithmetic Mean Calculation
load Tot_Rate_Var_rust.dat; % Obtain after Arithmetic Mean Calculation
load Tot_Rate_Var_selkov.dat; % Obtain after Arithmetic Mean Calculation
%% This calculation is for 100 points where total parameter variation between 0.01-0.015 range
%% For n=1 283K; n=2 288K; n=3 293K; n=4 298K; n=5 303K; n=6 308K; and n=7 313K;
n=7;
Luca = [Tot_Rate_Var_luca(n,:);Luca_osc_period(n,:)];
NFL = [Tot_Rate_Var_small_NFL(n,:);small_NFL_period(n,:)];
Selkov = [Tot_Rate_Var_selkov(n,:);subdep_period(n,:)];
Rust = [Tot_Rate_Var_rust(n,:);Rust_period(n,:)];
%%
A=Luca(Luca >=0.01 & Luca <= 0.015);
C = ismember(Luca,A);
columnsWithAllZeros = all(C == 0);
B = Luca(:, ~columnsWithAllZeros);
A1=NFL(NFL >=0.01 & NFL <= 0.015);
C1 = ismember(NFL,A1);
columnsWithAllZeros1 = all(C1 == 0);
B1 = NFL(:, ~columnsWithAllZeros1);
A2=Selkov(Selkov >=0.01 & Selkov <= 0.015);
C2 = ismember(Selkov,A2);
columnsWithAllZeros2 = all(C2 == 0);
B2 = Selkov(:, ~columnsWithAllZeros2);
A3=Rust(Rust >=0.01 & Rust <= 0.015);
C3 = ismember(Rust,A3);
columnsWithAllZeros3 = all(C3 == 0);
B3 = Rust(:, ~columnsWithAllZeros3);
%%
Bp=B(:,221:320);
B1p=B1(:,221:320);
B2p=B2(:,100:199);
B3p=B3(:,20:119);
figure(1)
hold on
plot(Bp(1,:),Bp(2,:),'o',...
    'color',[0, 0.4470, 0.7410],'MarkerSize',6,'DisplayName','Lucas Minimalistic Oscillatory Network')
plot(B1p(1,:),B1p(2,:),'or',...
    'MarkerSize',6,'DisplayName','Goodwins Negative Feedback Loop Model')
plot(B2p(1,:),B2p(2,:),'o',...
    'color',[0.9290, 0.6940, 0.1250],'MarkerSize',6,'DisplayName','Selkov Model')
plot(B3p(1,:),B3p(2,:),'o',...
    'color',[0.4660, 0.6740, 0.1880],'MarkerSize',6,'DisplayName','Rust Oscillatory Network')
xlabel('Total Parameter Variation') 
ylabel('Period of Oscillation') 
set(gca,'FontSize',10,'FontWeight','bold');
hAx=gca;
hAx.LineWidth=2; 
hLg=legend();
hLg.LineWidth=2;
%% save the files according to Temperature
% save -ascii luca_313K.dat Bp
% save -ascii NFL_313K.dat B1p
% save -ascii Selkov_313K.dat B2p
% save -ascii Rust_313K.dat B3p
%%
load luca_283K.dat;
load luca_288K.dat;
load luca_293K.dat;
load luca_298K.dat;
load luca_303K.dat;
load luca_308K.dat;
load luca_313K.dat;
Total_luca = [luca_283K(2,:);luca_288K(2,:);luca_293K(2,:);...
    luca_298K(2,:);luca_303K(2,:);luca_308K(2,:);luca_313K(2,:);];
% save -ascii data_luca.dat Total_luca;

load NFL_283K.dat;
load NFL_288K.dat;
load NFL_293K.dat;
load NFL_298K.dat;
load NFL_303K.dat;
load NFL_308K.dat;
load NFL_313K.dat;
Total_NFL = [NFL_283K(2,:);NFL_288K(2,:);NFL_293K(2,:);...
    NFL_298K(2,:);NFL_303K(2,:);NFL_308K(2,:);NFL_313K(2,:);];
% save -ascii data_NFL.dat Total_NFL;

load Selkov_283K.dat;
load Selkov_288K.dat;
load Selkov_293K.dat;
load Selkov_298K.dat;
load Selkov_303K.dat;
load Selkov_308K.dat;
load Selkov_313K.dat;
Total_Selkov = [Selkov_283K(2,:);Selkov_288K(2,:);Selkov_293K(2,:);...
    Selkov_298K(2,:);Selkov_303K(2,:);Selkov_308K(2,:);Selkov_313K(2,:);];
% save -ascii data_Selkov.dat Total_Selkov;

load Rust_283K.dat;
load Rust_288K.dat;
load Rust_293K.dat;
load Rust_298K.dat;
load Rust_303K.dat;
load Rust_308K.dat;
load Rust_313K.dat;
Total_Rust = [Rust_283K(2,:);Rust_288K(2,:);Rust_293K(2,:);...
    Rust_298K(2,:);Rust_303K(2,:);Rust_308K(2,:);Rust_313K(2,:);];
% save -ascii data_Rust.dat Total_Rust;
%% %CV calculations
load data_luca.dat;
load data_NFL.dat;
load data_Rust.dat;
load data_Selkov.dat;
CV_Luca=zeros(1,7);
CV_NFL=zeros(1,7);
CV_Selkov=zeros(1,7);
CV_Rust=zeros(1,7);
for l=1:7
CV_Luca(1,l)=100*(std(data_luca(l,:))/mean(data_luca(l,:)));
CV_NFL(1,l)=100*(std(data_NFL(l,:))/mean(data_NFL(l,:)));
CV_Selkov(1,l)=100*(std(data_Selkov(l,:))/mean(data_Selkov(l,:)));
CV_Rust(1,l)=100*(std(data_Rust(l,:))/mean(data_Rust(l,:)));
end
% save -ascii data_luca_cv.dat CV_Luca
% save -ascii data_NFL_cv.dat CV_NFL
% save -ascii data_Selkov_cv.dat CV_Selkov
% save -ascii data_Rust_cv.dat CV_Rust