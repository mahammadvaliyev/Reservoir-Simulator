%% Single-phase, 1D reservoir simulator

clear;
clc;
% specify dimensions
L=10500;   % ft
dy=400;    % ft
dz=40;     % ft
dx=500;    % ft
h=dz;      
ngrids=21;

% specify properties
phi=0.2;
k=15;                 % mD
mu=3;                 % cP
cfl=1.25*10^(-5);     % psi^-1, fluid compressibility
cform=0.75*10^(-5);   % psi^-1, formation compressibility
ct=cfl+cform;         % psi^-1
Bo=1.25;              % RB/STB
gamma=887.31*mu*Bo*dx/(k*dy*dz); % (cp*RB)/(mD*ft), useful const to define

% specify boundary conditions
P_i=3250;
% no flow BCs used in this project, so P0=P1, P22=P21

% well properties (same for all 3 wells)
rw=0.25;
re=0.14*(dx^2+dy^2)^(1/2);
J=7.08*10^(-3)*k*h/(mu*Bo*log(re/rw));  % bbl/psi, Productivity index


%% Leak-proof test

%specify simulation inputs
T=100; % days, simulation time
dt=10 ; % days, simulation timestep
ngrids=21;
beta=158.025*phi*mu*ct*(dx^2)/(k*dt)  % (cp*psi^-1*ft^2)/(mD*s), useful const to define
q4=0;
q11=0;
q21=0;

X_sol_leak = const_q(ngrids,T,dt,q4,q11,q21,P_i,gamma,beta);


figure(1)
plot(linspace(250,10250,21), X_sol_leak(:,1),'color','r','linewidth',2)
hold on
plot(linspace(250,10250,21), X_sol_leak(:,6),'color','b','linewidth',2)
plot(linspace(250,10250,21), X_sol_leak(:,11),'color','g','linewidth',2)
xlim([0 10500])
grid
legend('t=0', 't=50', 't=100', 'Location', 'NE')
title('Pressure distribution as function of x, t')
xlabel('Location (ft)') 
ylabel('Pressure (psi)') 
hold off

%% Symmetry test 

%specify simulation inputs
T=100; % days, simulation time
dt=10 ; % days, simulation timestep
ngrids=21;
beta=158.025*phi*mu*ct*(dx^2)/(k*dt)  % (cp*psi^-1*ft^2)/(mD*s), useful const to define
q4=0;
q11=100;
q21=0;

X_sol_sym = const_q(ngrids,T,dt,q4,q11,q21,P_i,gamma,beta);

figure(1)
plot(linspace(250,10250,21), X_sol_sym(:,1),'color','r','linewidth',2)
hold on
plot(linspace(250,10250,21), X_sol_sym(:,6),'-','color','b','linewidth',2)
plot(linspace(250,10250,21), X_sol_sym(:,11),'color','g','linewidth',2)
xline(5250,'--','color','k','linewidth',2)
grid
xlim([0 10500])
legend('t=0', 't=50', 't=100', 'Location', 'NE')
title('Pressure distribution as function of x, t')
xlabel('Location (ft)') 
ylabel('Pressure (psi)') 
hold off

%% Accuracy test

%specify simulation inputs
T=100; % days, simulation time
dt1=1 ; % days, simulation timestep
dt2=10 ; % days, simulation timestep
dt3=50 ; % days, simulation timestep
ngrids=21;

beta1=158.025*phi*mu*ct*(dx^2)/(k*dt1)
beta2=158.025*phi*mu*ct*(dx^2)/(k*dt2)
beta3=158.025*phi*mu*ct*(dx^2)/(k*dt3)

q4=0;
q11=100;
q21=0;

X_sol_acc1 = const_q(ngrids,T,dt1,q4,q11,q21,P_i,gamma,beta1);
X_sol_acc2 = const_q(ngrids,T,dt2,q4,q11,q21,P_i,gamma,beta2);
X_sol_acc3 = const_q(ngrids,T,dt3,q4,q11,q21,P_i,gamma,beta3);

figure(1)
plot(linspace(250,10250,21), X_sol_acc1(:,101),'color','r','linewidth',2)
hold on
plot(linspace(250,10250,21), X_sol_acc2(:,11),'-','color','b','linewidth',2)
plot(linspace(250,10250,21), X_sol_acc3(:,3),'color','g','linewidth',2)
grid
xlim([0 10500])
legend('dt=1', 'dt=10', 'dt=50', 'Location', 'NE')
title('Pressure distribution as function of x, t')
xlabel('Location (ft)') 
ylabel('Pressure (psi)') 
hold off

%% Non-uniform production rates in 3 wells

%specify simulation inputs
T=250; % days, simulation time
dt=10 ; % days, simulation timestep
ngrids=21;

beta=158.025*phi*mu*ct*(dx^2)/(k*dt)

q4=40;
q11=50;
q21=30;

X_sol_10a = const_q(ngrids,T,dt,q4,q11,q21,P_i,gamma,beta);

% 𝑝(𝑥,𝑡)vs. 𝑥for 𝑡=0,50,100,150,200,250 days
figure(1)
plot(linspace(250,10250,21), X_sol_10a(:,1),'color','r','linewidth',2)
hold on
plot(linspace(250,10250,21), X_sol_10a(:,6),'-','color','g','linewidth',2)
plot(linspace(250,10250,21), X_sol_10a(:,11),'color','b','linewidth',2)
plot(linspace(250,10250,21), X_sol_10a(:,16),'-','color','c','linewidth',2)
plot(linspace(250,10250,21), X_sol_10a(:,21),'color','y','linewidth',2)
plot(linspace(250,10250,21), X_sol_10a(:,26),'color','m','linewidth',2)
grid
xlim([0 10500])
legend('t=0', 't=50', 't=100','t=150','t=200','t=250', 'Location', 'SE', 'FontSize',8)
title('Pressure distribution as function of x, t')
xlabel('Location (ft)') 
ylabel('Pressure (psi)') 
hold off

%% 10. 𝑝𝑤𝑓(𝑡)vs time for all three wells in the same graph frame
figure(2)
plot(linspace(0,250,26), X_sol_10a(4,:)-q4/J,'color','g','linewidth',2)
hold on
plot(linspace(0,250,26), X_sol_10a(11,:)-q11/J,'color','y','linewidth',2)
plot(linspace(0,250,26), X_sol_10a(21,:)-q21/J,'color','m','linewidth',2)

grid
%xlim([0 10500])
legend('P_wf: Well Garnet','P_wf: Well Saph', 'P_wf: Well Opal', 'Location', 'NE', 'FontSize',8)
title('Pressure distribution as function of time')
xlabel('Time(days)') 
ylabel('Pressure (psi)') 
hold off

%% average pressure of the three simulation cells with wells.
figure(2)
plot(linspace(0,250,26), X_sol_10a(4,:),'color','r','linewidth',2)
hold on
plot(linspace(0,250,26), X_sol_10a(11,:),'color','b','linewidth',2)
plot(linspace(0,250,26), X_sol_10a(21,:),'color','c','linewidth',2)

grid
%xlim([0 10500])
legend('P_avg: Well Garnet', 'P_avg: Well Saph', 'P_avg: Well Opal','Location', 'NE', 'FontSize',8)
title('Pressure distribution as function of time')
xlabel('Time(days)') 
ylabel('Pressure (psi)') 
hold off

%% 𝑝𝑤f is constant and set to 1000 psia  for all three wells

%specify simulation inputs
T=250; % days, simulation time
dt=10 ; % days, simulation timestep
ngrids=21;

beta=158.025*phi*mu*ct*(dx^2)/(k*dt)

pwf4=1000;
pwf11=1000;
pwf21=1000;

X_sol = const_Pwf(ngrids,T,dt,pwf4,pwf11,pwf21,J,P_i,gamma,beta);

figure(2)
plot(linspace(0,250,26), (X_sol_10a(4,:)-pwf4)*J,'color','g','linewidth',2)
hold on
plot(linspace(0,250,26), (X_sol_10a(11,:)-pwf11)*J,'color','y','linewidth',2)
plot(linspace(0,250,26), (X_sol_10a(21,:)-pwf21)*J,'color','m','linewidth',2)

grid
%xlim([0 10500])
legend('q: Well Garnet','q: Well Saph', 'q: Well Opal', 'Location', 'NE', 'FontSize',8)
title('Flowrate as function time')
xlabel('Time(days)') 
ylabel('Flowrate (stb)') 
hold off

%% Cumulative flowrate as function time plot
figure(2)
plot(linspace(0,250,26), cumsum(X_sol_10a(4,:)),'color','r','linewidth',2)
hold on
plot(linspace(0,250,26),cumsum(X_sol_10a(11,:)),'color','b','linewidth',2)
plot(linspace(0,250,26),cumsum( X_sol_10a(21,:)),'color','c','linewidth',2)

grid
%xlim([0 10500])
legend('P_avg: Well Garnet', 'P_avg: Well Saph', 'P_avg: Well Opal','Location', 'SE', 'FontSize',8)
title('Cumulative flowrate as function time')
xlabel('Time(days)') 
ylabel('Cmulative flowrate (stb))') 
hold off
