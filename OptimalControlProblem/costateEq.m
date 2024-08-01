%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CoState equation 
%
% Author: Tuğba Akman Date: July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dp = costateEq(t,p,u1,Tu,x1,x2,x3,x4,x5,x6,x7,x8,x9,xt)

global Eq

dp = zeros(9,1);

Susc = interp1(xt,x1,t); 
Expo = interp1(xt,x2,t);
In = interp1(xt,x3,t); 
Is = interp1(xt,x4,t); 
Ih = interp1(xt,x5,t); 
Ihq = interp1(xt,x6,t); 
Q = interp1(xt,x7,t); 
Rec = interp1(xt,x8,t); 
Binc = interp1(xt,x9,t); 

u1 = interp1(Tu,u1,t); 

Npop = Susc  + Expo + In + Is + Ih + Ihq + Q + Rec;

R1S = -(Eq.betae*Expo + Eq.betae*In + Eq.betas*Is)*(((Npop-Ih-Q - Ihq-Susc)/((Npop-Ih-Ihq-Q)^2))) - Eq.q;
R1E = (Eq.betae*Expo + Eq.betae*In + Eq.betas*Is)*(((Susc)/((Npop-Ih-Ihq-Q)^2))) - Eq.betae*(Susc/(Npop-Ih-Ihq-Q));
R1In = (Eq.betae*Expo + Eq.betae*In + Eq.betas*Is)*(((Susc)/((Npop-Ih-Ihq-Q)^2))) - Eq.betae*(In/(Npop-Ih-Ihq-Q));
R1Is = (Eq.betae*Expo + Eq.betae*In + Eq.betas*Is)*(((Susc)/((Npop-Ih-Ihq-Q)^2))) - Eq.betas*(Is/(Npop-Ih-Ihq-Q));
R1Ih = 0;
R1Ihq = 0;
R1Q = 0;
R1R = (Eq.betae*Expo + Eq.betae*In + Eq.betas*Is)*(((Susc)/((Npop-Ih-Ihq-Q)^2))) ;
R1Bincr = 0;

R2S = (Eq.betae*Expo + Eq.betae*In + Eq.betas*Is)*(((Npop-Ih-Q - Ihq-Susc)/((Npop-Ih-Ihq-Q)^2)));
R2E = -R1E - Eq.k - Eq.q;
R2In = -R1In;
R2Is = -R1Is;
R2Ih = -R1Ih;
R2Ihq = 0;
R2Q = -R1Q;
R2R = -R1R;
R2Bincr = 0;

R3S = 0;
R3E = (1-Eq.rho)*Eq.k;
R3In = -Eq.gamman - Eq.q;
R3Is = 0;
R3Ih = 0;
R3Ihq = 0;
R3Q = 0;
R3R = 0;
R3Bincr = 0;

R4S = 0;
R4E = Eq.rho*Eq.k;
R4In = 0;
R4Is = -Eq.d - Eq.q - Eq.gammas; 
R4Ih = 0;
R4Ihq = 0;
R4Q = 0;
R4R = 0;
R4Bincr = 0;

R5S = 0;
R5E = 0;
R5In = 0;
R5Is = Eq.thetaa*(Eq.Bexist + Binc - Ih);
R5Ih = -Eq.gammah - Eq.nuh - Eq.thetaa*Is;
R5Ihq = 0;
R5Q = 0;
R5R = 0;
R5Bincr = Eq.thetaa*Is;

R6S = 0;
R6E = 0;
R6In = 0;
R6Is = Eq.d - Eq.thetaa*(Eq.Bexist + Binc - Ih);
R6Ih = Eq.thetaa*Is;
R6Ihq = -(Eq.gammahq + Eq.nuhq);
R6Q = 0;
R6R = 0;
R6Bincr = -Eq.thetaa*Is;

R7S = Eq.q;
R7E = Eq.q;
R7In = Eq.q;
R7Is = Eq.q;
R7Ih = 0;
R7Ihq = 0;
R7Q = -Eq.gammaq;
R7R = 0;
R7Bincr = 0;

R8S = 0;
R8E = 0;
R8In = Eq.gamman;
R8Is = Eq.gammas;
%R7Ih = Eq.phi;
R8Ih = Eq.gammah;
R8Ihq = Eq.gammahq;
R8Q = Eq.gammaq;
R8R = 0;
R8Bincr = 0;

R9S = 0;
R9E = 0;
R9In = 0;
R9Is = 0;
R9Ih = 0;
R9Ihq =(1+u1)*Eq.phi;
R9Q = 0;
R9R = 0;
R9Bincr = 0;

dp(1) = -p(1).*R1S - p(2).*R2S - p(3).*R3S - p(4).*R4S - p(5).*R5S - p(6).*R6S - p(7).*R7S - p(8).*R8S - p(9).*R9S;

dp(2) = -p(1).*R1E - p(2).*R2E - p(3).*R3E - p(4).*R4E - p(5).*R5E - p(6).*R6E - p(7).*R7E - p(8).*R8E  - p(9).*R9E;

dp(3) = -p(1).*R1In - p(2).*R2In - p(3).*R3In - p(4).*R4In - p(5).*R5In - p(6).*R6In - p(7).*R7In - p(8).*R8In - p(9).*R9In;
  
dp(4) = -p(1).*R1Is - p(2).*R2Is - p(3).*R3Is - p(4).*R4Is - p(5).*R5Is - p(6).*R6Is - p(7).*R7Is - p(8).*R8Is - p(9).*R9Is;
  
dp(5) = -p(1).*R1Ih - p(2).*R2Ih - p(3).*R3Ih - p(4).*R4Ih - p(5).*R5Ih - p(6).*R6Ih - p(7).*R7Ih - p(8).*R8Ih - p(9).*R9Ih;

dp(6) = -p(1).*R1Ihq - p(2).*R2Ihq - p(3).*R3Ihq - p(4).*R4Ihq - p(5).*R5Ihq - p(6).*R6Ihq - p(7).*R7Ihq - p(8).*R8Ihq - p(9).*R9Ihq -1;

dp(7) = -p(1).*R1Q - p(2).*R2Q - p(3).*R3Q - p(4).*R4Q - p(5).*R5Q - p(6).*R6Q - p(7).*R7Q - p(8).*R8Q - p(9).*R9Q ;

dp(8) = -p(1).*R1R - p(2).*R2R - p(3).*R3R - p(4).*R4R - p(5).*R5R - p(6).*R6R - p(7).*R7R - p(8).*R8R -  p(9).*R9R ;

dp(9) = -p(1).*R1Bincr - p(2).*R2Bincr - p(3).*R3Bincr - p(4).*R4Bincr - p(5).*R5Bincr - p(6).*R6Bincr - p(7).*R7Bincr - p(8).*R8Bincr - p(9).*R9Bincr;

end