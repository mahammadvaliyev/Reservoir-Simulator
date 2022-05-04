%% Two-phase, 1D reservoir simulator
clear;
clc;
%% Inputs

Pi=3000;
cf=1.25*10^(-5);
co=1.75*10^(-5);
cw=0.3*10^(-5);
k=100;
mu_o=10;
mu_w=1;

nx=100;
L=8000;
dz=50;
dy=200;
dx=L/nx;
h=dz; 
%%
rw=0.25;
s=0;
re=0.14*sqrt(dx^2+dy^2);

% Pressure dependent parameters
Soi=0.75; % function of P
Swc=0.25; % function of P
phi_i=0.2; % function of P
Boi=1.25;  % function of P
Bwi=1.05;   % function of P

% simulation timesteps
T=10000;
dt=10;
t_steps=T/dt; %n_time=T/dt+1; 

%incomplete constants

gamma_oil_wout_relp=887.3*(dx/(dy*dz))*(7.08/1000)*k*h/(mu_o*(log(re/rw)+s));
gamma_water_wout_relp=887.3*(dx/(dy*dz))*(7.08/1000)*k*h/(mu_w*(log(re/rw)+s));

%gamma_oil_wout_relp=0;
%gamma_water_wout_relp=0;
beta_wout_ct_phi=158*(dx^2)/dt; % %incomplete constant

%% Leak-proof test

% declare parameters and define their shapes
P_mat=zeros(nx,t_steps+1); % Pressure matrix
W=zeros(nx,t_steps);
C=zeros(nx,t_steps);
E=zeros(nx,t_steps);
RHS=zeros(nx,t_steps);

Bo_avg=zeros(nx,t_steps+1);
Bo_i=zeros(nx,t_steps+1);
Bw_avg=zeros(nx,t_steps+1);
Bw_i=zeros(nx,t_steps+1);


phi=zeros(nx,t_steps+1);
Swup=zeros(nx,t_steps+1);
krw=zeros(nx,t_steps+1);
kro=zeros(nx,t_steps+1);
mob_w=zeros(nx,t_steps+1);
mob_o=zeros(nx,t_steps+1);
ct=zeros(nx,t_steps+1);

beta=zeros(nx,t_steps+1);
gamma_w=zeros(1,t_steps+1);
gamma_o=zeros(1,t_steps+1);

So=zeros(nx,t_steps+1);
Sw=zeros(nx,t_steps+1);

%%
% initialize parameters
P_mat(:,1)=ones(nx,1)*Pi;
So(:,1)=ones(nx,1)*Soi;
Sw(:,1)=ones(nx,1)*Swc;
ct(:,1)=cf*ones(nx,1)+co*So(:,1)+cw*Sw(:,1);

Bo_avg(:,1)=Boi*ones(nx,1);
Bo_i(:,1)=Boi*ones(nx,1);
Bw_avg(:,1)=Bwi*ones(nx,1);
Bw_i(:,1)=Bwi*ones(nx,1);

phi(:,1)=phi_i*ones(nx,1);
Swup(:,1)=Swc*ones(nx,1);

krw(:,1)=ones(nx,1)*rel_perm_water(Swc);
kro(:,1)=ones(nx,1)*rel_perm_oil(Swc);
mob_w(:,1)=Bw_i(:,1).*k.*krw(:,1)./(mu_w.*Bw_avg(:,1));
mob_o(:,1)=Bo_i(:,1).*k.*kro(:,1)./(mu_o.*Bo_avg(:,1));

beta(:,1)=beta_wout_ct_phi.*phi(:,1).*ct(:,1);
gamma_o(1)=gamma_oil_wout_relp*kro(20,1);
gamma_w(1)=gamma_water_wout_relp*krw(20,1);
%gamma_o(1)=0;
%gamma_w(1)=0;
Pwf=1000;

AQ=1;
qo=zeros(1,t_steps+1);
qw=zeros(1,t_steps+1);
%% Solve for P and S
% Use old Bo, Bw, S, krw, mob to solve for P
% After updating P, update: phi, Bo, Bw, So, Sw, krw, mob, ct, beta, gamma
for i=2:t_steps+1
    for j=1:nx
        if (AQ==1) 
            Swup(1,i-1)=1; % or 0.25
            krw(1,i-1)=rel_perm_water(Swup(1,i-1));
            kro(1,i-1)=rel_perm_oil(Swup(1,i-1));
            mob_w(1,i-1)=Bw_i(1,i-1)*k*krw(1,i-1)/(mu_w*Bw_avg(1,i-1));
            mob_o(1,i-1)=Bo_i(1,i-1)*k*kro(1,i-1)/(mu_o*Bo_avg(1,i-1));
        else
            Swup(1,i-1)=0.25;
            krw(1,i-1)=0;
            kro(1,i-1)=0;
            mob_w(1,i-1)=0;
            mob_o(1,i-1)=0;
        end  
        if (j==1)
            W(j,i-1)=mob_w(j,i-1)+mob_o(j,i-1);
            E(j,i-1)=mob_w(j+1,i-1)+mob_o(j+1,i-1);
            C(j,i-1)=-(W(j,i-1)+E(j,i-1)+beta(j,i-1));
            if (AQ==1)
                RHS(j,i-1)=-beta(j,i-1)*P_mat(j,i-1)-W(j,i-1)*Pi; %%%%
            else
                RHS(j,i-1)=-beta(j,i-1)*P_mat(j,i-1)-W(j,i-1)*P_mat(j,i-1); %%%%
            end
        elseif (j==100)
            W(j,i-1)=mob_w(j,i-1)+mob_o(j,i-1);
            E(j,i-1)=0;
            C(j,i-1)=-(W(j,i-1)+E(j,i-1)+beta(j,i-1)+gamma_o(i-1)+gamma_w(i-1));
            RHS(j,i-1)=-beta(j,i-1)*P_mat(j,i-1)-(gamma_o(i-1)+gamma_w(i-1))*Pwf;
        else
            W(j,i-1)=mob_w(j,i-1)+mob_o(j,i-1);
            E(j,i-1)=mob_w(j+1,i-1)+mob_o(j+1,i-1);
            C(j,i-1)=-(W(j,i-1)+E(j,i-1)+beta(j,i-1));
            RHS(j,i-1)=-beta(j,i-1)*P_mat(j,i-1);
        end
    end
    P_mat(:,i)=tridiag(C(:,i-1), W(:,i-1), E(:,i-1), RHS(:,i-1));
    phi(:,i)=phi(:,i-1).*exp(cf*ones(nx,1).*(P_mat(:,i)-P_mat(:,i-1)));
    for j=1:nx
        if (j==1)
            if (AQ==1)
                Bo_avg(j,i)=Boi*exp(-co*(0.5*(Pi+P_mat(j,i))-Pi));
                Bw_avg(j,i)=Bwi*exp(-cw*(0.5*(Pi+P_mat(j,i))-Pi));
                Bo_i(j,i)=Boi;
                Bw_i(j,i)=Bwi;
                So(j,i)=(Bo_i(j,i)/phi(j,i))* ( (phi(j,i-1)*So(j,i-1)/Bo_i(j,i-1)) + (mob_o(j,i-1)*Pi-(mob_o(j,i-1)+mob_o(j+1,i-1))*P_mat(j,i)+mob_o(j+1,i-1)*P_mat(j+1,i))/(158*Bo_i(j,i-1)*dx^2/dt));
                Sw(j,i)=(Bw_i(j,i)/phi(j,i))* ( (phi(j,i-1)*Sw(j,i-1)/Bw_i(j,i-1)) + (mob_w(j,i-1)*Pi-(mob_w(j,i-1)+mob_w(j+1,i-1))*P_mat(j,i)+mob_w(j+1,i-1)*P_mat(j+1,i))/(158*Bw_i(j,i-1)*dx^2/dt));
                Swup(1,i-1)=1; % or 0.25           
            else
                Bo_avg(j,i)=Boi*exp(-co*(P_mat(j,i)-Pi));
                Bw_avg(j,i)=Bwi*exp(-cw*(P_mat(j,i)-Pi));
                Bo_i(j,i)=Bo_i(j,i-1)*exp(-co*(P_mat(j,i)-P_mat(j,i-1)));
                Bw_i(j,i)=Bw_i(j,i-1)*exp(-cw*(P_mat(j,i)-P_mat(j,i-1)));
                So(j,i)=(Bo_i(j,i)/phi(j,i)) * ( (phi(j,i-1)*So(j,i-1)/Bo_i(j,i-1)) + (-mob_o(j+1,i-1)*P_mat(j,i)+mob_o(j+1,i-1)*P_mat(j+1,i))/(158*Bo_i(j,i-1)*dx^2/dt));
                Sw(j,i)=(Bw_i(j,i)/phi(j,i)) * ( (phi(j,i-1)*Sw(j,i-1)/Bw_i(j,i-1)) + (-mob_w(j+1,i-1)*P_mat(j,i)+mob_w(j+1,i-1)*P_mat(j+1,i))/(158*Bw_i(j,i-1)*dx^2/dt));
                Swup(j,i)=sat_up(P_mat(j,i),P_mat(j,i),Sw(j,i),Sw(j,i));
                krw(j,i)=rel_perm_water(Swup(j,i));
                kro(j,i)=rel_perm_oil(Swup(j,i));
                mob_w(j,i)=Bw_i(j,i)*k*krw(j,i)/(mu_w*Bw_avg(j,i));
                mob_o(j,i)=Bo_i(j,i)*k*kro(j,i)/(mu_o*Bo_avg(j,i));
            end
            %Swup(j,i)=sat_up(P_mat(j,i),P_mat(j,i),Sw(j,i),Sw(j,i));
            %krw(j,i)=rel_perm_water(Swup(j,i));
            %kro(j,i)=rel_perm_oil(Swup(j,i));
            %mob_w(j,i)=Bw_i(j,i)*k*krw(j,i)/(mu_w*Bw_avg(j,i));
            %mob_o(j,i)=Bo_i(j,i)*k*kro(j,i)/(mu_o*Bo_avg(j,i));
        elseif (j==100)
            Bo_avg(j,i)=Boi*exp(-co*(0.5*(P_mat(j,i)+P_mat(j-1,i))-Pi));
            Bw_avg(j,i)=Bwi*exp(-cw*(0.5*(P_mat(j,i)+P_mat(j-1,i))-Pi));
            Bo_i(j,i)=Bo_i(j,i-1)*exp(-co*(P_mat(j,i)-P_mat(j,i-1)));
            Bw_i(j,i)=Bw_i(j,i-1)*exp(-cw*(P_mat(j,i)-P_mat(j,i-1)));             
            So(j,i)=(Bo_i(j,i)/phi(j,i))* ( (phi(j,i-1)*So(j,i-1)/Bo_i(j,i-1)) + (mob_o(j,i-1)*P_mat(j-1,i)-(mob_o(j,i-1))*P_mat(j,i)-gamma_o(i-1)*(P_mat(j,i)-Pwf))/(158*Bo_i(j,i-1)*dx^2/dt));
            Sw(j,i)=(Bw_i(j,i)/phi(j,i))* ( (phi(j,i-1)*Sw(j,i-1)/Bw_i(j,i-1)) + (mob_w(j,i-1)*P_mat(j-1,i)-(mob_w(j,i-1))*P_mat(j,i)-gamma_w(i-1)*(P_mat(j,i)-Pwf))/(158*Bw_i(j,i-1)*dx^2/dt));
            Swup(j,i)=sat_up(P_mat(j-1,i),P_mat(j,i),Sw(j-1,i),Sw(j,i));
            krw(j,i)=rel_perm_water(Swup(j,i));
            kro(j,i)=rel_perm_oil(Swup(j,i));
            mob_w(j,i)=Bw_i(j,i)*k*krw(j,i)/(mu_w*Bw_avg(j,i));
            mob_o(j,i)=Bo_i(j,i)*k*kro(j,i)/(mu_o*Bo_avg(j,i));
        else
            Bo_avg(j,i)=Boi*exp(-co*(0.5*(P_mat(j,i)+P_mat(j-1,i))-Pi));
            Bw_avg(j,i)=Bwi*exp(-cw*(0.5*(P_mat(j,i)+P_mat(j-1,i))-Pi));
            Bo_i(j,i)=Bo_i(j,i-1)*exp(-co*(P_mat(j,i)-P_mat(j,i-1)));
            Bw_i(j,i)=Bw_i(j,i-1)*exp(-cw*(P_mat(j,i)-P_mat(j,i-1)));
            So(j,i)=(Bo_i(j,i)/phi(j,i))* ( (phi(j,i-1)*So(j,i-1)/Bo_i(j,i-1)) + (mob_o(j,i-1)*P_mat(j-1,i)-(mob_o(j,i-1)+mob_o(j+1,i-1))*P_mat(j,i)+mob_o(j+1,i-1)*P_mat(j+1,i))/(158*Bo_i(j,i-1)*dx^2/dt));
            Sw(j,i)=(Bw_i(j,i)/phi(j,i))* ( (phi(j,i-1)*Sw(j,i-1)/Bw_i(j,i-1)) + (mob_w(j,i-1)*P_mat(j-1,i)-(mob_w(j,i-1)+mob_w(j+1,i-1))*P_mat(j,i)+mob_w(j+1,i-1)*P_mat(j+1,i))/(158*Bw_i(j,i-1)*dx^2/dt));
            Swup(j,i)=sat_up(P_mat(j-1,i),P_mat(j,i),Sw(j-1,i),Sw(j,i));
            krw(j,i)=rel_perm_water(Swup(j,i));
            kro(j,i)=rel_perm_oil(Swup(j,i));
            mob_w(j,i)=Bw_i(j,i)*k*krw(j,i)/(mu_w*Bw_avg(j,i));
            mob_o(j,i)=Bo_i(j,i)*k*kro(j,i)/(mu_o*Bo_avg(j,i));
        end
        disp(So(j,i)+Sw(j,i))
    end  
    ct(:,i)=cf*ones(nx,1)+co*So(:,i)+cw*Sw(:,i);
    beta(:,i)=beta_wout_ct_phi.*phi(:,i).*ct(:,i);
    gamma_o(:,i)=gamma_oil_wout_relp*kro(20,i);
    gamma_w(:,i)=gamma_water_wout_relp*krw(2,i);
    qo(:,i)=(P_mat(20,i)-Pwf)*gamma_o(:,i)/Bo_avg(20,i);
    qw(:,i)=(P_mat(20,i)-Pwf)*gamma_w(:,i)/Bw_avg(20,i);
end           
   
%% Pressures
figure(1)
plot(linspace(0,8000,100), P_mat(:,1),'color','r','linewidth',2)
hold on
plot(linspace(0,8000,100), P_mat(:,6),'color','b','linewidth',2)
plot(linspace(0,8000,100), P_mat(:,11),'color','g','linewidth',2)
xlim([-100 8100])
grid
legend('t=0', 't=50', 't=100', 'Location', 'NE')
title('Pressure distribution as function of x, t')
xlabel('Location (ft)') 
ylabel('Pressure (psi)') 
hold off
%% Flowrates
figure(2)
plot(linspace(0,T,t_steps), qo(2:t_steps+1),'color','r','linewidth',2)
hold on
plot(linspace(0,T,t_steps), qw(2:t_steps+1),'color','b','linewidth',2)
xlim([-100 10100])
grid
legend('oil rate', 'water rate', 'Location', 'NE')
title('Oil & water Cmulative flowrates, stb/d')
xlabel('Time (days)') 
ylabel('Flowrate (stb/day)') 
hold off

%% Cumulative flowrate as function time plot
figure(2)
plot(linspace(0,T,t_steps), cumsum(qo(2:t_steps+1)),'color','r','linewidth',2)
hold on
plot(linspace(0,T,t_steps),cumsum(qw(2:t_steps+1)),'color','b','linewidth',2)
grid
xlim([-100 10100])
legend('oil rate', 'water rate', 'Location', 'SE')
title('Oil & water Cmulative flowrates, stb/d')
xlabel('Time (days)') 
ylabel('Flowrate (stb/day)') 
hold off

%% Water cut
figure(3)
plot(linspace(0,T,t_steps), qw(2:t_steps+1)./(qw(2:t_steps+1)+qo(2:t_steps+1)),'color','r','linewidth',2)

grid
xlim([-100 10100])
legend('well water cut', 'Location', 'SE')
title('Water cut vs time ')
xlabel('Time (days)') 
ylabel('water cut ') 

%% Water saturation
N=L*dy*dz*phi_i*Soi;
Np=cumsum(qo(2:t_steps+1));
RF=round(Np/N,2);
RF_vec=0:0.05:0.45;
indeces=zeros(1,1);
for i=0:0.05:0.45
    c=find(RF==i,1);
    indeces = [indeces c];
end
%%
figure(4)
colors=["red","green","blue","cyan","magenta","yellow","black","b"]
for i=1:length(indeces)
    plot(linspace(0,8000,100), Sw(:,indeces(i)+1),'color',colors(i),'linewidth',2)
    ylim([0 1])
    hold on
    title('Water saturation distribution as function of x, t')
    xlabel('Location (ft)') 
    ylabel('Water saturation') 
    legend('RF=0', 'RF=0.05','RF=0.1','RF=0.15','RF=0.2','RF=0.25','RF=0.3','RF=0.35','Location', 'NE')
    %hold off
end
%%

title('Water saturation distribution as function of x, t')
xlabel('Location (ft)') 
ylabel('Water saturation (psi)') 
hold off


