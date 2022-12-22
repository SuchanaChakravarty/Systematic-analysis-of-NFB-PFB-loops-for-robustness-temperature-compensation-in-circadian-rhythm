% ************************************************************************
% Gillespie algorithim for Selkov-like positive feedback oscillator with 
% an additional negative feedback loop model
% ************************************************************************
% Rate constants declare
% ************************************************************************
clear
clc
t = 0.0; % initial time
t_end = 1150.0; % final time
t_sample = 0.02; % sampling time
Nr = 100; % number of trajectory
Vs = 100; % scaling factor
Temp = 298; % temperature in Kelvin
t_run = t_end/t_sample+1;

% ************************************************************************
% Set-1 : Small k3 Large k5 Case
% ************************************************************************
k1 = 383.83.*exp(-(16.35.*10^3)./(8.3144598.*Temp));
k2 = 383.83.*exp(-(27.5.*10^3)./(8.3144598.*Temp));
k3 = 383.83.*exp(-(12.0203.*10^3)./(8.3144598.*Temp));
k4 = 383.83.*exp(-(14.8691.*10^3)./(8.3144598.*Temp));
k5 = 383.83.*exp(-(5.17.*10^3)./(8.3144598.*Temp));

% ************************************************************************
% Set-2 : Large k3 small k5 Case
% ************************************************************************
% k1 = 383.83.*exp(-(21.5.*10^3)/(8.3144598.*Temp));
% k2 = 383.83.*exp(-(27.5.*10^3)./(8.3144598.*Temp));
% k3 = 383.83.*exp(-(5.17*10^3)/(8.3144598.*Temp));
% k4 = 383.83.*exp(-(14.8691.*10^3)./(8.3144598.*Temp));
% k5 = 383.83.*exp(-(12.0203*10^3)/(8.3144598.*Temp));

% ************************************************************************
% Set-3 : Large k3 Large k5 Case
% ************************************************************************
% k1 = 383.83.*exp(-(20.877.*10^3)/(8.3144598.*Temp));
% k2 = 383.83.*exp(-(27.5.*10^3)./(8.3144598.*Temp));
% k3 = 383.83.*exp(-(5.17*10^3)/(8.3144598.*Temp));
% k4 = 383.83.*exp(-(14.8691.*10^3)./(8.3144598.*Temp));
% k5 = 383.83.*exp(-(5.17*10^3)/(8.3144598.*Temp));

% ************************************************************************
% Initialization
% ************************************************************************
j = 1.0; % counter for time 
nY = 0.48.*Vs; % initial Y represented in numbers
nU = 0.0075.*Vs; % initial U represented in numbers

% ************************************************************************
% loop for trajectory starts
% ************************************************************************
for q=1:Nr    
% ************************************************************************
%  Initial states before starting SSA
% ************************************************************************
q
t_array(j,q)=t;  
Y_array(j,q)=nY; 
U_array(j,q)=nU;
% ************************************************************************
% Reaction Proepensities
% ************************************************************************
while t < t_end
a1 = ((k1.*Vs.^2)./(Vs+(k5.*nU))); % in numbers
a2 = k2.*nY + k3.*((nU.^2)./(Vs.^2)).*nY; % in numbers
a3 = k2.*nY + k3.*((nU.^2)./(Vs.^2)).*nY; % in numbers
a4 = k4.*nU; % in numbers
h = [a1 a2 a3 a4];
%combined rxn propensities
    h0 = sum(h);
    if h0<=0
        h0=0.000000001;
    end
    r1=0+rand*(1-0);
            if (r1<=0)
                r1=0.00000001;
            end
    %time update
     t_next = ((1/h0)*(log(1/r1)));
     t = t + t_next;
     %determine next reaction
       i=1; mu=0; amu=0; r2=rand;
         while amu < r2*h0
              mu = mu + 1;
              amu = amu + h(i); 
              i = i + 1;
         end
%reactions
    if mu == 1
        nY = nY + 1;
    elseif mu == 2
        nY = nY - 1;
        elseif mu == 3
        nU = nU + 1;

    elseif mu == 4
        nU = nU - 1;
    end
% avoid neg or zero values of the species    
if nY<=0
nY=0.000001;
end
if nU<=0
nU=0.000001;
end
%store/output time and species 
if t >= j*t_sample
     j=j+1;
     t_array(j,q)=t;
     Y_array(j,q)=nY;
     U_array(j,q)=nU; 
end
end
trow{:,q}=t_array(:,q); % store time in a cell as the size differes in each run
yi{:,q}=U_array(:,q); % store species in cell as the legth differs in each run
% smooth the U_array in order to avoid small peaks while calculating period of oscillation
for ik=1:15
    yk{q,1}=yi{:,q};
    ym{q,ik}=smooth(yk{q,ik});
    yk{q,ik+1}=ym{q,ik};
end
% period of oscillation calculations
w=(t_run-1)/t_end;
[peakval,locval]=findpeaks(ym{q,15},'MinPeakHeight',0.2*Vs,'MinPeakDistance',0.8*w);
% [peakval,locval]=findpeaks(ym{q,15},'MinPeakHeight',1*Vs,'MinPeakDistance',0.8*w); % for lk3-sk5, vs=100
pkval{:,q}=peakval; % storing peakval in cell
lcval{:,q}=locval; % storing locval in cell
tval{:,q}=trow{1,q}(lcval{:,q}); % identifying the exact time from lcval and stored in a cell
period(:,q) = mean(diff(tval{:,q})); % period calculation
% ************************************************************************
% SSA ends here
% ************************************************************************
j=1.0; % reset time counter
t=0.0; % reset initial time
nY=Y_array(1,q); % nY update
nU=U_array(1,q); % nU update
end    
% ************************************************************************
% trajectory loop ends
% ************************************************************************
%%
s=std(period);
m=mean(period);
%%
figure(1)
hold on;
plot(t_array,U_array)
%%
figure(2)
plot(period)
%%
U_array_smooth=ym;
std_period=s;
mean_period=m;
% save sk3_lk5_Vs10.mat U_array U_array_smooth t_array trow std_period mean_period period
CV_model=100*(std_period/mean_period);
