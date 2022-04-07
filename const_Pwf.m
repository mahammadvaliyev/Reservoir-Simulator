function [X_sol] = const_Pwf(ngrids,T,dt,pwf4,pwf11,pwf21,J,P_i,gamma,beta)

% Function returns Pressure matrix (ngrids x ntime) that contains pressure
% values at each grid cells (rows) for each timestep (columns)

% Inputs:
% ngrids - number of grid cells
% T - simulation time
% dt - timestep size
% q4, q11, q21 - flowrate values for wells at cells 4, 11, 21
% beta - useful constant that simplifies the code

% Output:
% X_sol - pressure values matrix

n_time=T/dt+1;  

X_sol=zeros(ngrids,n_time);
X_sol(:,1)=P_i;

%populate coefficients: a, b, c
b=ones(ngrids,1);
c=ones(ngrids,1);
a=zeros(ngrids,1);

for i=1:ngrids
    if (i==1 || i==21)
        a(i,1)=-(1+beta);
    elseif (i==4 || i==11 || i==20)
        a(i,1)=-(2+beta+J*gamma);
    else
        a(i,1)=-(2+beta);
    end
end

% vector that will have non-zero values at cells with wells
well_term1=gamma*pwf4*J;
well_term2=gamma*pwf11*J;
well_term3=gamma*pwf21*J;
q_RHS=zeros(ngrids,1);

for i=1:ngrids
    if (i==4)
        q_RHS(i,1)=well_term1;
    elseif (i==11)
        q_RHS(i,1)=well_term2;
    elseif (i==20)
        q_RHS(i,1)=well_term3;
    else
        q_RHS(i,1)=0;
    end
end

f=zeros(ngrids,1);

for i=1:n_time-1
    
    f(:,1)=-beta*X_sol(:,i)+q_RHS;
    X_sol(:,i+1)=tridiag(a,b,c,f);
    
end


end

