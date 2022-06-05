load subdep_period.dat; % get that value after the simulation from Arithmatic-Mean folder
load subdep_period_avg.dat; % get that value after the simulation from Time-Course folder
A = subdep_period;
B = subdep_period_avg;
C = zeros(7,100);
C(1,1:100) = A(1,451:550) - B(1,1);
C(2,1:100) = A(2,451:550) - B(1,2);
C(3,1:100) = A(3,451:550) - B(1,3);
C(4,1:100) = A(4,451:550) - B(1,4);
C(5,1:100) = A(5,451:550) - B(1,5);
C(6,1:100) = A(6,451:550) - B(1,6);
C(7,1:100) = A(7,451:550) - B(1,7);
D=C.^2;
SSE1=sum(D,2);
n=100; % number of observation
k=4; % number of parameter
BICk=n*log(SSE1)-n*log(n)+k*log(n);
Temp = [283; 288; 293; 298; 303; 308; 313];
data = [Temp BICk];
% save -ascii bic_subdep.dat data
%%