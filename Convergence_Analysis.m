clc; clear; close all;

threshold = 1e-6;
iterations = 200;
L = 1;

auc(iterations) = 0;
for i=1:iterations
C_i = Solve(i);
dx_i = L/(length(C_i)-1);
% auc(i) = Simpson3_8(C_i, dx_i);
auc(i) = trapz(C_i) *dx_i;
end


%Relative to N-1 Elements -> Meaningless
% error(1,iterations-1) = 0; %Pre-allocates memory for 99 floating-pt #. Equivalent to zeros(1,iteraitons-1)
% for i=1:iterations-1
%    error(i) = abs((auc(i)-auc(i+1)))/auc(i);
% end


%Relative to Analytical solution
x = linspace(0,L,10000);
C_actual = (1/(1-exp(5)))*(1-exp(5*x));
auc_actual = integral(@(x)(1/(1-exp(5)))*(1-exp(5*x)),0,1);
error_ac = abs(auc-auc_actual)/auc_actual;


%Plot
f = figure( 'NumberTitle','off',...
            'Name','Convergence Analysis | Joe Issa',...
            'Menubar','figure');
f.Position(3) = 700;

% subplot(1,2,1)
% loglog(2:iterations,error,'b')
% title(["Error relative to the previous"; "solution for N-1 elements"])
% xlabel("Number of Elements {\itN}")
% ylabel("Relative error")
% hold on
% grid minor
% plot([1,iterations],[threshold threshold],'r')
% N_conv = find(error<threshold,1)+1;
% if ~isempty(N_conv)
% plot([N_conv,N_conv],[1e-10,error(N_conv-1)],'m--')
% else
% disp('Maybe increase the value of \iterations\?')    
% end
% xlim([0 iterations])


% subplot(1,2,2)
loglog(1:iterations,error_ac,'b')
title("Error relative to the analytical solution")
xlabel("Number of Elements {\itN}")
ylabel("Relative error")
hold on
grid minor
plot([1,iterations],[threshold threshold],'r')
N_conv_ac = find(error_ac<threshold,1);
if ~isempty(N_conv_ac)
plot([N_conv_ac,N_conv_ac],[1e-8,error_ac(N_conv_ac)],'m--')
else
disp('Maybe increase the value of \iterations\?')    
end
% plot(N_conv_ac, error_ac(N_conv_ac),'m','MarkerSize',30)
xlim([0 iterations])



% function auc = Simpson3_8(C_int, dx)
% auc = 3*dx/8*( C_int(1)+C_int(end)+3*sum(C_int(2:3:end-2))+...
%            3*sum(C_int(3:3:end-1))+2*sum(C_int(4:3:end-3)) );
% end

function C_inter = Solve(in)
D = 2; %m^2/s, Diffusivity of the Material
u = 10; %m/s, Flow Velocity
L = 1; %m, Length

Elts = in;
order = 3;

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
       row = B(:,i);
       col = B(:,j);
       Ki = sparse(row,col,K_e(i,j),Nodes,Nodes,N_nz);
       K = K + Ki;
    end
end

Q(1:Nodes,1) = [linspace(0,0,Elts)';1;linspace(0,0,Nodes-Elts-1)'];
K(1,:) = [1,linspace(0,0,Nodes-1)];
K(Elts+1,:) = [linspace(0,0,Elts)';1;linspace(0,0,Nodes-Elts-1)'];

C = K\Q;

%POST-PROCESSING
%===============

%Interpolate inter elements
%
% N = cast(cast(250/Elts,"int64")+1,'like',1); % number of internal nodes evaluated and plotted, per Element
if Elts<1000
N = 50;
C_inter(1,(N+1)*Elts+1) = 0;
x_inter(1, (N+1)*Elts+1) = 0;

inter(order+1,N) = 0;
step = 2/(N+1);

for i=1:order+1
    inter(i,:) = polyval(psi(i,:),-1+step:step:1-step);
end

i=1;
y = [C(B(i,1)),C(B(i,:))'*inter,C(B(i,end))];
x = linspace((i-1)*h_e,i*h_e,N+2);
C_inter(1,1:N+2) = y;
x_inter(1,1:N+2) = x;
for i=2:Elts
% display(N+3+(i-2)*(N+1):N+2+(i-1)*(N+1))
y = [C(B(i,1)),C(B(i,:))'*inter,C(B(i,end))];
x = linspace((i-1)*h_e,i*h_e,N+2);
C_inter(1,N+3+(i-2)*(N+1):N+2+(i-1)*(N+1)) = y(2:end);
x_inter(1,N+3+(i-2)*(N+1):N+2+(i-1)*(N+1)) = x(2:end);
end
else
    C_inter = C(Mesh);
end

end