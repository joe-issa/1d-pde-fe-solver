function [C, B, psi, Mesh] = SolveForC (param)

D = param.D; %m^2/s, Diffusivity of the Material
u = param.u; %m/s, Flow Velocity
L = param.L; %m, Length

Elts = param.Elts;
order = param.order;
Nodes = Elts * order +1; % = Nodes + (order-1) * Nodes + 1
h_e = L/Elts;
dzeta_over_dx = 2/h_e;

Mesh(1:order:Nodes) = 1:Elts+1;
for i=2:order
Mesh(i:order:Nodes-order+i-1) = Elts+i:order-1:Nodes-order+i;
end

psi = importdata('psi_coef.dat');
dpsi = importdata('dpsi_coef.dat');

psi = cell2mat(psi(order));
dpsi = cell2mat(dpsi(order));

%Gauss Points
GP.Abs = importdata('GP_Abs.dat');
GP.Wei = importdata('GP_Wei.dat');

N_GP = order;
GP.Abs = cell2mat(GP.Abs(N_GP));
GP.Wei = cell2mat(GP.Wei(N_GP));

psi_val(order+1, N_GP) = 0;
dpsi_val(order+1, N_GP) = 0;


%PRE-PROCESSING
%==============
Elts = Elts * ((Elts>0)-(Elts<0)); % -> Absolute Value
order = (mod(cast(order,"int64")-1,5)+1); %order can only be from +1 to +5


%RESOLUTION
%==========
K_e(order,order) = 0;

for i=1:order+1
    psi_val(i,:) = polyval(psi(i,:),GP.Abs);
    dpsi_val(i,:) = polyval(dpsi(i,:),GP.Abs);
end

for i=1:order+1
    for j=1:order+1
    K_e(i,j) = (dzeta_over_dx*D*dpsi_val(j,:).*dpsi_val(i,:)+...
        u*psi_val(i,:).*dpsi_val(j,:))*GP.Wei;
    end
end


B(Elts,order+1) = 0;

B(:,1) = 1:Elts;
B(:,order+1) = 2:Elts+1;
for i=2:order
    B(:,i) = Elts+i:order-1:Nodes-order+i;
end


N_nz = Elts*((order+1)^2 -1) + 1; %Before imposing any boundary condition
K = spalloc(Nodes,Nodes,N_nz);
for i=1:order+1
    for j=1:order+1
%        row = B(:,i);
%        col = B(:,j);
%        Ki = sparse(row,col,K_e(i,j),Nodes,Nodes,N_nz);
       K = K + sparse(B(:,i),B(:,j),K_e(i,j),Nodes,Nodes,N_nz);
    end
end


Q(1:Nodes,1) = [linspace(0,0,Elts)';1;linspace(0,0,Nodes-Elts-1)'];
K(1,:) = [1,linspace(0,0,Nodes-1)];
K(Elts+1,:) = [linspace(0,0,Elts)';1;linspace(0,0,Nodes-Elts-1)'];

C = K\Q;
end