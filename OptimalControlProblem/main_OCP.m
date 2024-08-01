%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to simulate optimal intervention strategies
%
% Author: Tuğba Akman Date: July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main_OCP()

close all
clc
format long

set(0, 'defaultaxesfontsize',16)

global Eq

%Initial conditions
init_Susc = 1e+03;
init_Expo = 1e+03; 
init_In = 2*1e+04;  
init_Is = 1.5*1e+04;  
init_Ih = 1e+03;
init_Ihq = 0;
init_Q = 1e+05;
init_Rec = 1e+05;
init_Binc = 0; 

% Parameters
Eq.Bexist= init_Ih;
Eq.k = 1/14;  
Eq.gammas = 1/2.9;  
Eq.gamman = 1/7; 
Eq.betae = +7.8e-05;
Eq.betas = 9.8e+02*0.01;
Eq.gammah = 0.36; 
Eq.nuh = 0.0068;
Eq.nuq = 0.0034;
Eq.phi = 0.1; 
Eq.rho = 0.0078;
Eq.thetaa = 2.4e-02;
Eq.q = 0.05;
Eq.q = 0.05;
Eq.q = 0.05;
Eq.gammaq = 0.18;
Eq.d = 1;
Eq.phi0 = 0.0001;
Eq.gammahq = 0.0018; 
Eq.nuhq = 0.068;

% Weights in the cost functional
Eq.weight1 = 1;  

% Control constraints
M1=0;    % Lower bound for control
M2=1; % Upper bound for control

% Time discretization
t0=0;                   % initial time
tf = 25;                % final time
Interval=200;           % Number of subintervals
dt = (tf-t0)/Interval; % Temporal increment
Tu=t0:(1*dt):tf;        % ınterpolation points

% day to start intervention
intervention_starts = 0; % days

%ICs
initx=[init_Susc, init_Expo, init_In, init_Is, init_Ih,init_Ihq, init_Q, init_Rec, init_Binc]'; % IC for state
initlambda=[0,0,0,0,0,0,0,0,0]';                  % FC for adjoint

%Allocate for uncontrolled case - MUST BE EQUAL TO ZERO!!!!
u1=0*ones(size(Tu))';

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

%No control
[Tx_uncon,X_uncon] = ode15s(@(t,x) stateEq_no_control(t,x,u1,Tu), [t0, tf], ...
    initx, options);

%%
% Time discretization for intervention
t0_treatment=intervention_starts;           % treatment starts at t=t0_treatment
Interval_treatment=200;                 % Number of subintervals
dt_treatment = ( tf-t0_treatment)/Interval_treatment; % Temporal increment
Tu_treatment=t0_treatment:(1*dt_treatment):tf;        % interpolation points

% Allocate for intervention
u1=0*ones(size(Tu_treatment))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solve OCP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initx_treatment = initx;

theta_convex=linspace(0.01,0.99,100);

lambda1=zeros(size(Tu_treatment));
lambda2=zeros(size(Tu_treatment));
lambda3=zeros(size(Tu_treatment));
lambda4=zeros(size(Tu_treatment));
lambda5=zeros(size(Tu_treatment));
lambda6=zeros(size(Tu_treatment));
lambda7=zeros(size(Tu_treatment));
lambda8=zeros(size(Tu_treatment));
lambda9=zeros(size(Tu_treatment));
x1=zeros(size(Tu_treatment));
x2=zeros(size(Tu_treatment));
x3=zeros(size(Tu_treatment));
x4=zeros(size(Tu_treatment));
x5=zeros(size(Tu_treatment));
x6=zeros(size(Tu_treatment));
x7=zeros(size(Tu_treatment));
x8=zeros(size(Tu_treatment));
x9=zeros(size(Tu_treatment));

nn=1; test = 2;
while(test > 1e-5)

    
    oldu1 = u1;
    oldx1=x1;
    oldx2=x2;
    oldx3=x3;
    oldx4=x4;
    oldx5=x5;
    oldx6=x6;
    oldx7=x7;
    oldx8=x8;
    oldx9=x9;
    oldlambda1=lambda1;
    oldlambda2=lambda2;
    oldlambda3=lambda3;
    oldlambda4=lambda4;
    oldlambda5=lambda5;
    oldlambda6=lambda6;
    oldlambda7=lambda7;
    oldlambda8=lambda8;
    oldlambda9=lambda9;

    % solve the state eqn forward
    [Tx,X] = ode15s(@(t,x) stateEq_control(t,x,u1,Tu_treatment), Tu_treatment, ...
        initx_treatment, options);

    
    % solve the adjoint eqn backward
    x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);  x5 = X(:,5);  x6 = X(:,6);  x7 = X(:,7); x8 = X(:,8); x9 = X(:,9);
        [Tlambda,Lambda] = ode15s(@(t,lambda) costateEq(t,lambda,u1,Tu_treatment,x1,x2,x3,x4,x5,x6,x7,x8,x9,Tx), ...
        flip(Tu_treatment), initlambda, options);

    lambda1 = Lambda(:,1); lambda2 = Lambda(:,2); lambda3 = Lambda(:,3);
    lambda4 = Lambda(:,4); lambda5 = Lambda(:,5); lambda6 = Lambda(:,6); 
    lambda7 = Lambda(:,7); lambda8 = Lambda(:,8); lambda9 = Lambda(:,9);

    % Interpolation
    lambda1 = interp1(Tlambda,lambda1,Tx);
    lambda2 = interp1(Tlambda,lambda2,Tx);
    lambda3 = interp1(Tlambda,lambda3,Tx);
    lambda4 = interp1(Tlambda,lambda4,Tx);
    lambda5 = interp1(Tlambda,lambda5,Tx);
    lambda6 = interp1(Tlambda,lambda6,Tx);
    lambda7 = interp1(Tlambda,lambda7,Tx);
    lambda8 = interp1(Tlambda,lambda8,Tx);
    lambda9 = interp1(Tlambda,lambda9,Tx);

    % Optimality condition
    Cond = Eq.phi*x6;
    optCond1 = -lambda9.*Cond;  

    % Project the control
    u1 = min(M2, max(M1, optCond1));
   
    % Update the control
    for dd=1:1:length(theta_convex)

        thetaa=theta_convex(dd);
        u1_convex_for_inside(:,dd) = (thetaa)*u1 + (1-thetaa)*oldu1;
        J_for_inside(dd)=cost(Tx,X,u1_convex_for_inside(:,dd));
        
    end

    [~, pos_min_for_inside] = min(J_for_inside);

    u1_convex_for_inside(:,pos_min_for_inside);
    u1 = u1_convex_for_inside(:,pos_min_for_inside);
    
    %Use same interpolation  points
    Tu_treatment=Tx;

    % Stopping criteria on relative error
    temp1 = norm((oldu1 - u1))/norm((u1));
    
    temp2 = norm((oldx1 - x1))/norm(x1);
    temp3 = norm((oldx2 - x2))/norm(x2);
    temp4 = norm((oldx3 - x3))/norm(x3);
    temp5 = norm((oldx4 - x4))/norm(x4);
    temp6 = norm((oldx5 - x5))/norm(x5);
    temp7 = norm((oldx6 - x6))/norm(x6);
    temp8 = norm((oldx7 - x7))/norm(x7);
    temp17 = norm((oldx8 - x8))/norm(x8);
    temp18 = norm((oldx9 - x9))/norm(x9);    

    temp9 = norm((oldlambda1 - lambda1))/norm((lambda1));
    temp10 = norm((oldlambda2 - lambda2))/norm((lambda2));
    temp11 = norm((oldlambda3 - lambda3))/norm((lambda3));
    temp12 = norm((oldlambda4 - lambda4))/norm((lambda4));
    temp13 = norm((oldlambda5 - lambda5))/norm((lambda5));
    temp14 = norm((oldlambda6 - lambda6))/norm((lambda6));
    temp15 = norm((oldlambda7 - lambda7))/norm((lambda7));
    temp16 = norm((oldlambda8 - lambda8))/norm((lambda8));
    temp19 = norm((oldlambda9 - lambda9))/norm((lambda9));

    test = max(temp1, max(temp2, max( temp3, max( temp4, max(temp5,...
        max(temp6, max(temp7, max(temp8, max(temp9, max( temp10, max(temp11,...
        max(temp12, max(temp13, max(temp14, max(temp15, max(temp16, max(temp17, max(temp18,temp19))))))))))))))))))

    nn=nn+1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% solve the state eqn forward
[Tx,X] = ode15s(@(t,x) stateEq_control(t,x,u1,Tu_treatment), Tu_treatment, ...
        initx_treatment, options);


figure(1)
plot(Tx,M2*ones(size(Tx)),'-k','LineWidth',2)
hold on
plot(Tx,u1,'--r','LineWidth',2)
title('Control function')
xlabel('t')
ylabel('u(t)')
legend('Maximum intervention','u(t)','Location','SouthWest')
%xlim([10,tf])
%xline(t0_treatment,'--k','HandleVisibility','off')
grid on

figure(211)
plot(Tx,X(:,1), '-r','LineWidth',2);
hold on
plot(Tx_uncon,X_uncon(:,1), '--b','LineWidth',2);
title('S(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
legend('With control','No control','Location','NorthEast')
grid on
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Susc_ocp','eps')
saveas(fig,'Susc_ocp','fig')

figure(212)
plot(Tx,X(:,2), '-r','LineWidth',2);
hold on
plot(Tx_uncon,X_uncon(:,2), '--b','LineWidth',2);
title('E(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
legend('With control','No control','Location','NorthEast')
grid on
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Expo_ocp','eps')
saveas(fig,'Expo_ocp','fig')

figure(213)
plot(Tx,X(:,3), '-r','LineWidth',2);
hold on
plot(Tx_uncon,X_uncon(:,3), '--b','LineWidth',2);
title('I_n(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
legend('With control','No control','Location','NorthEast')
grid on
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'In_ocp','eps')
saveas(fig,'In_ocp','fig')

figure(214)
plot(Tx,X(:,4), '-r','LineWidth',2);
hold on
plot(Tx_uncon,X_uncon(:,4), '--b','LineWidth',2);
title('I_s(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
legend('With control','No control','Location','NorthEast')
grid on
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Is_ocp','eps')
saveas(fig,'Is_ocp','fig')

figure(215)
plot(Tx,X(:,5), '-r','LineWidth',2);
hold on
plot(Tx_uncon,X_uncon(:,5), '--b','LineWidth',2);
title('I_h(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
legend('With control','No control','Location','NorthEast')
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Ih_ocp','eps')
saveas(fig,'Ih_ocp','fig')

figure(216)
plot(Tx,X(:,6), '-r','LineWidth',2);
hold on
plot(Tx_uncon,X_uncon(:,6), '--b','LineWidth',2);
title('I_{hq}(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
legend('With control','No control','Location','NorthEast')
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Ihq_ocp','eps')
saveas(fig,'Ihq_ocp','fig')


figure(217)
plot(Tx,X(:,7), '-r','LineWidth',2);
hold on
plot(Tx_uncon,X_uncon(:,7), '--b','LineWidth',2);
title('Q(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
legend('With control','No control','Location','NorthEast')
grid on
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Q_ocp','eps')
saveas(fig,'Q_ocp','fig')

figure(218)
plot(Tx,X(:,8), '-r','LineWidth',2);
hold on
plot(Tx_uncon,X_uncon(:,8), '--b','LineWidth',2);
title('R(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
legend('With control','No control','Location','SouthEast')
grid on
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Rec_ocp','eps')
saveas(fig,'Rec_ocp','fig')


figure(219)
plot(Tx,X(:,9), '-r','LineWidth',2);
hold on
plot(Tx_uncon,X_uncon(:,9), '--b','LineWidth',2);
title('B_{incr}(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
legend('With control','No control','Location','SouthEast')
grid on
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Bincr_ocp','eps')
saveas(fig,'Bincr_ocp','fig')

figure(220)
plot(Tx,Eq.rho*Eq.k*X(:,2), '-r','LineWidth',2);
hold on
plot(Tx_uncon,Eq.rho*Eq.k*X_uncon(:,2), '--b','LineWidth',2);
legend('With control','No control','Location','SouthWest')
title('Incidences')
xlabel('Time (days)','fontweight','normal','fontsize',18)
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Incidences_ocp','eps')
saveas(fig,'Incidences_ocp','fig')

figure(221)
plot(Tx,Eq.nuh*X(:,5)+Eq.nuhq*X(:,6), '-r','LineWidth',2);
hold on
plot(Tx_uncon,Eq.nuh*X_uncon(:,5)+Eq.nuhq*X_uncon(:,6), '--b','LineWidth',2);
legend('With control','No control','Location','NorthEast')
title('Deaths')
xlabel('Time (days)','fontweight','normal','fontsize',18)
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Deaths_ocp','eps')
saveas(fig,'Deaths_ocp','fig')

end

%Sub-function

function J=cost(Tx,X,u1)
global Eq

J=(trapz(Tx,X(:,6)) + Eq.weight1*0.5*trapz(Tx,u1.^2));
end

