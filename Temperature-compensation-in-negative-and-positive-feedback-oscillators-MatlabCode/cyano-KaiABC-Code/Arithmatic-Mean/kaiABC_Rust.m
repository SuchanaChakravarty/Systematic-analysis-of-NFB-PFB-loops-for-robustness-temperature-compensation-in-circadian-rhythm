function dcdt = kaiABC_Rust(t,c,Temp,...
    Ae_khalf,Ae_kUTA,Ae_kDTA,...
    Ae_kTU0,Ae_kTUA,Ae_kTDA,Ae_kSDA,...
    Ae_kDS0,Ae_kDSA,Ae_kUSA,Ae_kSU0,...
    Ae_kSUA)
% variables
T = c(1);
ST = c(2);
S = c(3);
Tot_Rate_Var = c(4);
% KaiA = 1.3;
KaiC = 3.4;
% parameters, random pre-exponential factor(Ae)
% Generate Random Number first before Running this code
khalf = Ae_khalf*exp(-(16.8330449*10^3)/(8.3144598*Temp));
kUT0 = 0.0; 
kUTA = Ae_kUTA*exp(-(16.56527851*10^3)/(8.3144598*Temp));
kDT0 = 0.0; 
kDTA = Ae_kDTA*exp(-(19.08885839*10^3)/(8.3144598*Temp));
kTU0 = Ae_kTU0*exp(-(18.60866545*10^3)/(8.3144598*Temp));
kTUA = Ae_kTUA*exp(-(21.0046209*10^3)/(8.3144598*Temp));
kTD0 = 0.0; 
kTDA = Ae_kTDA*exp(-(18.5744178*10^3)/(8.3144598*Temp));
kSD0 = 0.0; 
kSDA = Ae_kSDA*exp(-(16.43132499*10^3)/(8.3144598*Temp));
kDS0 = Ae_kDS0*exp(-(17.64373845*10^3)/(8.3144598*Temp));
kDSA = Ae_kDSA*(-exp(-(17.56984493*10^3)/(8.3144598*Temp)));
kUS0 = 0.0; 
kUSA = Ae_kUSA*exp(-(22.00905957*10^3)/(8.3144598*Temp));
kSU0 = Ae_kSU0*exp(-(20.21073081*10^3)/(8.3144598*Temp));
kSUA = Ae_kSUA*(-exp(-(19.73888331*10^3)/(8.3144598*Temp)));

U = max( 0 ,KaiC - T - ST - S);
A = max( 0 ,1.3 - 2.*S);

% Parameters defined as function
kUT = kUT0 + ((kUTA.*A)./(khalf+A));
kDT = kDT0 + ((kDTA.*A)./(khalf+A)); 
kTU = kTU0 + ((kTUA.*A)/(khalf+A));
kTD = kTD0 + ((kTDA.*A)./(khalf+A)); 
kSD = kSD0 + ((kSDA.*A)./(khalf+A)); 
kDS = kDS0 + ((kDSA.*A)./(khalf+A));
kUS = kUS0 + ((kUSA.*A)./(khalf+A));
kSU = kSU0 + ((kSUA.*A)./(khalf+A));
% parameter as base value at each temperature, pre-exponential factor
% =Ae=383.83
khalfp = 383.83*exp(-(16.8330449*10^3)/(8.3144598*Temp)); 
kUTAp = 383.83*exp(-(16.56527851*10^3)/(8.3144598*Temp));
kDTAp = 383.83*exp(-(19.08885839*10^3)/(8.3144598*Temp));
kTU0p = 383.83*exp(-(18.60866545*10^3)/(8.3144598*Temp));
kTUAp = 383.83*exp(-(21.0046209*10^3)/(8.3144598*Temp));
kTDAp = 383.83*exp(-(18.5744178*10^3)/(8.3144598*Temp));
kSDAp = 383.83*exp(-(16.43132499*10^3)/(8.3144598*Temp));
kDS0p = 383.83*exp(-(17.64373845*10^3)/(8.3144598*Temp));
kDSAp = 383.83*(-exp(-(17.56984493*10^3)/(8.3144598*Temp))); 
kUSAp = 383.83*exp(-(22.00905957*10^3)/(8.3144598*Temp));
kSU0p = 383.83*exp(-(20.21073081*10^3)/(8.3144598*Temp));
kSUAp = 383.83*(-exp(-(19.73888331*10^3)/(8.3144598*Temp)));

kUTp = kUT0 + ((kUTAp.*A)./(khalfp+A));
kDTp = kDT0 + ((kDTAp.*A)./(khalfp+A)); 
kTUp = kTU0p + ((kTUAp.*A)./(khalfp+A));
kTDp = kTD0 + ((kTDAp.*A)./(khalfp+A)); 
kSDp = kSD0 + ((kSDAp.*A)./(khalfp+A)); 
kDSp = kDS0p + ((kDSAp.*A)./(khalfp+A));
kUSp = kUS0 + ((kUSAp.*A)./(khalfp+A));
kSUp = kSU0p + ((kSUAp.*A)./(khalfp+A));
%total no of parameter
n=12;
%total parameter variation
Z=((abs(khalf-khalfp)/khalfp)+(abs(kUTA-kUTAp)/kUTAp)+(abs(kDTA-kDTAp)/kDTAp)...
    +(abs(kTU0-kTU0p)/kTU0p)+(abs(kTUA-kTUAp)/kTUAp)+(abs(kTDA-kTDAp)/kTDAp)...
    +(abs(kSDA-kSDAp)/kSDAp)+(abs(kDS0-kDS0p)/kDS0p)+(abs(kDSA-kDSAp)/kDSAp)...
    +(abs(kUSA-kUSAp)/kUSAp)+(abs(kSU0-kSU0p)/kSU0p)+(abs(kSUA-kSUAp)/kSUAp))/n;
% initialization 
dcdt = zeros(4,1);
% differential equations
dcdt(1) = (kUT.*U + kDT.*ST - kTU.*T - kTD.*T);
dcdt(2) = (kTD.*T + kSD.*ST - kDT.*ST - kDS.*ST);
dcdt(3) = (kUS.*U + kDS.*ST - kSU.*S - kSD.*S);
dcdt(4) = (Z- Tot_Rate_Var);