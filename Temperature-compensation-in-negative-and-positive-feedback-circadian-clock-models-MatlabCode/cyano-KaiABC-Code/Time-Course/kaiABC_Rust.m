function dcdt = kaiABC_Rust(t,c,Temp)
% variables
T = c(1);
ST = c(2);
S = c(3);
% KaiA = 1.3;
KaiC = 3.4;
% parameters
khalf = 383.83*exp(-(16.8330449*10^3)/(8.3144598*Temp));
kUT0 = 0.0; 
kUTA = 383.83*exp(-(16.56527851*10^3)/(8.3144598*Temp));
kDT0 = 0.0; 
kDTA = 383.83*exp(-(19.08885839*10^3)/(8.3144598*Temp));
kTU0 = 383.83*exp(-(18.60866545*10^3)/(8.3144598*Temp));
kTUA = 383.83*exp(-(21.0046209*10^3)/(8.3144598*Temp));
kTD0 = 0.0; 
kTDA = 383.83*exp(-(18.5744178*10^3)/(8.3144598*Temp));
kSD0 = 0.0; 
kSDA = 383.83*exp(-(16.43132499*10^3)/(8.3144598*Temp));
kDS0 = 383.83*exp(-(17.64373845*10^3)/(8.3144598*Temp));
kDSA = 383.83*(-exp(-(17.56984493*10^3)/(8.3144598*Temp)));
kUS0 = 0.0; 
kUSA = 383.83*exp(-(22.00905957*10^3)/(8.3144598*Temp));
kSU0 = 383.83*exp(-(20.21073081*10^3)/(8.3144598*Temp));
kSUA = 383.83*(-exp(-(19.73888331*10^3)/(8.3144598*Temp)));

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
% initialization 
dcdt = zeros(3,1);
% differential equations
dcdt(1) = (kUT.*U + kDT.*ST - kTU.*T - kTD.*T);
dcdt(2) = (kTD.*T + kSD.*ST - kDT.*ST - kDS.*ST);
dcdt(3) = (kUS.*U + kDS.*ST - kSU.*S - kSD.*S);