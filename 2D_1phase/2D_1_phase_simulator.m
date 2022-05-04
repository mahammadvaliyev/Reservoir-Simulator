%% Define inputs

clear;
clc;

Nx=11;
Ny=9;
dlt_x=500;
dlt_y=500;
dlt_z=100;
dlt_t=10;

c_o=1.75e-05;
c_f=1.25e-05;
c_t=c_f+c_o;
Bo=1.25;

m=Nx*Ny;

lambda_S=zeros(m,1);
lambda_N=zeros(m,1);
lambda_W=zeros(m,1);
lambda_E=zeros(m,1);
lambda_C=zeros(m,1);

rhv=zeros(m,1);

p_ini=3500;

pwf=1500;

phi_ini=0.2;
miu=5;
kx=10;
ky=2;

re=0.28*sqrt(ky*dlt_x^2+kx*dlt_y^2)/(sqrt(kx)+sqrt(ky));
rw=0.25;
S=1;

for i=1:m
pv(i)=3500;
end

%% Initial Pressure distribution
pv=pv.';
p2d=reshape(pv,Nx,Ny)';
surf(p2d);
szTitle=sprintf('Reservoir Pressure Distribution t=0 days');
axis([0 Nx 0 Ny 2000 3600])
title(szTitle,'fontsize',15,'fontweight','bold','Color','blue');
xlabel('x-direction','fontsize',10,'fontweight','bold')
ylabel('y-direction','fontsize',10,'fontweight','bold')
%% Solve for Pressure distribution at timesteps 1...50

P_matrix=zeros(m,50);
for j=1:50
    Jw=0.00708*sqrt(kx*ky)*dlt_z/(miu*Bo*exp(c_o*(p_ini-pv(50)))*(log(re/rw)+S));
    for i=1:m
    alpha(i)=158.03*c_t*phi_ini*exp(c_f*(pv(i)-p_ini));
    end
    for i=1:m
    lambda_S(i)=ky/(miu);
    lambda_N(i)=ky/(miu);
    lambda_W(i)=kx/(miu);
    lambda_E(i)=kx/(miu);
    end
    for i=1:Nx
    lambda_S(i)=0;
    lambda_N((Ny-1)*Nx+i)=0;
    end
    for i=1:Ny
    lambda_W(1+(i-1)*Nx)=0;
    lambda_E(i*Nx)=0;
    end
    for i=1:m
    lambda_C(i)=-lambda_E(i)-lambda_W(i)-lambda_N(i)-lambda_S(i)-(alpha(i)*dlt_x^2/dlt_t);
    end
    Bin=[lambda_S lambda_W lambda_C lambda_E lambda_N];
    d=[-Nx -1 0 1 Nx];
    sparse = spdiags(Bin,d,m,m+1);
    matrix=full(sparse);
    matrix(:,m+1)=[];
    for i=1:m
    rhv(i)=-alpha(i)*pv(i)*dlt_x^2/dlt_t;
    end
    q(j)=(pv(50)-pwf)*Jw;
    rhv(50)=-alpha(50)*pv(50)*dlt_x^2/dlt_t+(887.3*dlt_x*q(j)*Bo*exp(c_o*(p_ini-pv(50)))/(dlt_z*dlt_y));
    pv_new=matrix\rhv;
    P_matrix(:,j)=pv_new;
    pv=pv_new;
    
    
end
%% Plot P(x,t)
m=50 %specify timestep
p2d=reshape(P_matrix(:,m),Nx,Ny)';
surf(p2d);
szTitle=sprintf('Reservoir Pressure Distribution t=%.0f days',m*dlt_t);
axis([0 Nx 0 Ny 2000 3600])
%title(szTitle,'fontsize',15,'fontweight','bold');
title(szTitle,'fontsize',15, 'Color','blue');
xlabel('x-direction','fontsize',10,'fontweight','bold')
ylabel('y-direction','fontsize',10,'fontweight','bold')
     

%% Plot Q
q(51)=(pv(50)-pwf)*Jw;
for k=1:51
    time_v(k)=(k-1)*dlt_t;
end
figure(2);
plot(time_v,q,'Color','red');
xlabel('time,days');
ylabel('q,STB/D');
title('Well Production Rate','Color','blue');