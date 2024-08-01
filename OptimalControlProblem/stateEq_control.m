%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State equation 
%
% Author: Tuğba Akman Date: July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = stateEq_control(t,y,u1,Tu)

global Eq

dy = zeros(9,1);

u1 = interp1(Tu,u1,t); 

Susc = y(1);
Expo = y(2);
In = y(3);
Is = y(4);
Ih = y(5);
Ihq = y(6);
Q = y(7);
Rec = y(8);
Binc = y(9);

Npop = Susc  + Expo + In + Is + Ih + Ihq + Q + Rec;

dy(1) = -((Eq.betae*Expo*Susc + Eq.betae*In*Susc + Eq.betas*Is*Susc)/(Npop-Ih-Q - Ihq)) - Eq.q*Susc;
dy(2) = ((Eq.betae*Expo*Susc + Eq.betae*In*Susc + Eq.betas*Is*Susc)/(Npop-Ih-Q - Ihq)) - Eq.k*Expo - Eq.q*Expo;
dy(3) = (1-Eq.rho)*Eq.k*Expo - Eq.gamman*In - Eq.q*In;
dy(4) = Eq.rho*Eq.k*Expo - (Eq.d + Eq.q + Eq.gammas)*Is;
dy(5) = Eq.thetaa*(Eq.Bexist + Binc - Ih)*Is - (Eq.gammah + Eq.nuh)*Ih;
dy(6) = (Eq.d - Eq.thetaa*(Eq.Bexist + Binc - Ih))*Is - (Eq.gammahq + Eq.nuhq)*Ihq;
dy(7) = Eq.q*Susc + Eq.q*Expo + Eq.q*In + Eq.q*Is - Eq.gammaq*Q;
dy(8) = Eq.gamman*In + Eq.gammas*Is + Eq.gammah*Ih + Eq.gammahq*Ihq + Eq.gammaq*Q;
dy(9) = (1+u1)*Eq.phi*Ihq; 
end
