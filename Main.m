%Close figures, uifigure fui, clear workspace and command window, stop
%previous diary
close all;
if exist('fui','var') == 1
    delete(fui)
    diary off
end
clear; clc;

%PHYSICAL PARAMETERS: Transport of Chemical Species
%===================
D = 2; %m^2/s, Diffusivity of the Material
u = 10; %m/s, Flow Velocity
L = 2; %m, Length


%MESH
%====

%INPUT PARAMETERS
order = 3;
want_to_plot = true; %Interactive GUI for Elts <= 1e3
want_to_log = true; %Takes more time to run

Elts = input('# of Elements (must be numerical, 29 taken by default): '); %Keep below 1e3 for interactive gui
if isempty(Elts) || Elts==0
   Elts = 29; % From Convergence Analysis with Analytical Solution
end


%============================================================
%============================================================

if want_to_log==1
    start_date = datestr(now,'dd-mm-yy-HHMMSS');
    Log(start_date) 
end

%LogToConsole
disp("Solving for the concentrations, with parameters:")
fprintf(1,"\t"+Elts+" Elements\n")
fprintf(1,"\tMaterial Diffusivity of "+D+" m^2/s\n")
fprintf(1,"\tFlow Velocity of "+u+" m/s\n")
fprintf(1,"\tDistance between the parallel plates of "+L+" meter(s)\n")
disp("Using:")
fprintf(1,"\tStandard Galerkin Formulation\n")
fprintf(1,"\tSolution approximated using Lagrange shape functions of order "+order+"\n")

tic
Nodes = Elts * order +1; % = Nodes + (order-1) * Nodes + 1
h_e = L/Elts;
dzeta_over_dx = 2/h_e; %Ratio of legnths of elements

Mesh(1:order:Nodes) = 1:Elts+1; %Fills the all physial nodes & allocates all needed memory
% For virtual nodes:
for i=2:order
Mesh(i:order:Nodes-order+i-1) = Elts+i:order-1:Nodes-order+i;
end

% %Calculating psi and dpsi for any polynomial order

% zetta = -1:2/order:1;
% n = nchoosek(zetta,order);
% 
% psi(order+1,order+1) = 0; 
% dpsi(order+1,order) = 0; 
% 
% for i=1:order+1
%     p = poly(n(i,:));
%     psi(i,:) = p/polyval(p,zetta(end+1-i));
%     dpsi(i,:) = polyder(psi(i,:));
% end
% psi = psi(end:-1:1,:);
% dpsi = dpsi(end:-1:1,:);

%For a faster code, we will use stored values of psi and dpsi
psi = importdata('psi_coef.dat');
dpsi = importdata('dpsi_coef.dat');

psi = cell2mat(psi(order));
dpsi = cell2mat(dpsi(order));

%Gauss Points
%A struct of arrays was tested to be faster than 2x1-D arrays or 1x2-D
%array in the Gausss Quadrature process. 
%Coefficients from: https://pomax.github.io/bezierinfo/legendre-gauss.html
GP.Abs = importdata('GP_Abs.dat');
GP.Wei = importdata('GP_Wei.dat');

%In this application, N_GP always simplifies to 'order'
N_GP = order; %((order+order-1)+1)/2 = 2(order)/2 = order
GP.Abs = cell2mat(GP.Abs(N_GP));
GP.Wei = cell2mat(GP.Wei(N_GP));

psi_val = zeros(order+1, N_GP);
dpsi_val = zeros(order+1, N_GP);


%PRE-PROCESSING (Quick fixes to avoid errors)
%==============
Elts =  abs(Elts); %OR Elts * ((Elts>0)-(Elts<0));  -> Absolute Value
order = (mod(cast(order,"int64")-1,5)+1); %order can only be from +1 to +5
% clear p n zetta; %Free up memory space to store more values

%RESOLUTION
%==========
K_e = zeros(order+1,order+1);

% 1. Evaluate shape functions at all Gauss Points.
for i=1:order+1
    psi_val(i,:) = polyval(psi(i,:),GP.Abs);
    dpsi_val(i,:) = polyval(dpsi(i,:),GP.Abs);
end

% 2. Find [K^e]. In our case, it's equal for all elements so we calculate
% only once
for i=1:order+1
    for j=1:order+1
    K_e(i,j) = (dzeta_over_dx*D*dpsi_val(j,:).*dpsi_val(i,:)+...
        u*psi_val(i,:).*dpsi_val(j,:))*GP.Wei;
    end
end

% 3. Derive the Connectivity Matrix B
% B(1:Elts,1:order+1)=[1:Elts ; Elts+2:order-1:Nodes-1 ; Elts+3:order-1:Nodes ; 2:Elts+1]';
% One line method, slow and specific to order = 3

%The method below is ~faster~ for high numbers of Elements (Elts)
%where speed is a concern and ~general~ for any order
B = zeros(Elts,order+1);

B(:,1) = 1:Elts;
B(:,order+1) = 2:Elts+1;
for i=2:order
    B(:,i) = Elts+i:order-1:Nodes-order+i;
end

% 4. Assemble!
N_nz = Elts*((order+1)^2 -1) + 1; %Before imposing any boundary condition
K = spalloc(Nodes,Nodes,N_nz);
for i=1:order+1
    for j=1:order+1
%        row = B(:,i);
%        col = B(:,j);
%        Ki = sparse(B(:,i),B(:,j),K_e(i,j),Nodes,Nodes,N_nz);
       K = K + sparse(B(:,i),B(:,j),K_e(i,j),Nodes,Nodes,N_nz);
    end
end

% 5. Impose Boundary Conditions
Q(1:Nodes,1) = [linspace(0,0,Elts)';1;linspace(0,0,Nodes-Elts-1)'];
K(1,:) = [1,linspace(0,0,Nodes-1)];
K(Elts+1,:) = [linspace(0,0,Elts)';1;linspace(0,0,Nodes-Elts-1)'];

% 6. Solve matrix system and Done!
C = K\Q;
disp("Resolution Complete!")
toc
if want_to_log==1 %Save matrix
    cd Logs
    writematrix([linspace(0,1,Nodes)',C], start_date+"_results.txt") 
    cd ..
end




%POST-PROCESSING
%===============
if want_to_plot==1
screen_size = get(0,'ScreenSize');
screen_size = screen_size(3:4);

%Interactive GUI for a small numbers of Elements (less than 1K)
% if Elts < 10001 && want_to_plot == true
disp("Building GUI... First time takes more time")
tic
fui = uifigure('Name','Convection-Diffusion 1D FE Solver | Joe Issa');

%Containers: panels and grids
g1 = uigridlayout(fui,[2,2]); %Divide into 4

p_par = uipanel(g1,'Title','Parameters');
p_par.Layout.Row = 1;
p_par.Layout.Column = 1;

g3 = uigridlayout(p_par,[10,4]);

p_hm = uipanel(g1,'Title','Heat Map');
p_hm.Layout.Row = 2;
p_hm.Layout.Column = 1;

%=========================================

%Plot Area
C_plot = uiaxes(g1);
C_plot.Layout.Row = [1,2];
C_plot.Layout.Column = 2;

grid(C_plot,'on')
C_plot.Title.String = ["Concentration of cations as a function";"of the position between the plates"];
C_plot.XLabel.String = "Distance from 0 to {\itL} (in m)";
C_plot.YLabel.String = "Cation Concentration (in mol/m^3)";
%=========================================

%Parameters Area
 %Material Diffusivity
  %Edit Field
DD = uieditfield(g3,'numeric');
DD.Value = D;
DD.Limits = [1e-15, 1e5];
DD.Layout.Row = 1;
DD.Layout.Column = 2;
  %Label
Lbl_DD = uilabel(g3,'Text','Diffusivity (m2/s): ');
Lbl_DD.HorizontalAlignment = 'right';
Lbl_DD.FontSize = 14;
Lbl_DD.Layout.Row = 1;
Lbl_DD.Layout.Column = 1;
 %Length L
  %Edit Field
LL = uieditfield(g3,'numeric');
LL.Value = L;
LL.Limits = [1e-2, 1e5];
LL.Layout.Row = 2;
LL.Layout.Column = 2;
  %Label
Lbl_LL = uilabel(g3,'Text','Length (m): ');
Lbl_LL.HorizontalAlignment = 'right';
Lbl_LL.FontSize = 14;
Lbl_LL.Layout.Row = 2;
Lbl_LL.Layout.Column = 1;

 %Flow Velocity
  %Edit Field
uu = uieditfield(g3,'numeric');
uu.Value = u;
uu.Limits = [-1e5, 1e5];
uu.Layout.Row = 1;
uu.Layout.Column = 4;
  %Label
Lbl_uu = uilabel(g3,'Text','Velocity (m/s): ');
Lbl_uu.HorizontalAlignment = 'right';
Lbl_uu.FontSize = 14;
Lbl_uu.Layout.Row = 1;
Lbl_uu.Layout.Column = 3;

 %Enter Number of Elements
  %Edit Field
N_Elts = uieditfield(g3,'numeric');
N_Elts.Value = Elts;
N_Elts.Limits = [1, 1e6]; %Limited to 1K to avoid losing GUI responsiveness
N_Elts.Layout.Row = 4;
N_Elts.Layout.Column = 2;
  %Label
Lbl_N_Elts = uilabel(g3,'Text','# of Elements: ');
Lbl_N_Elts.HorizontalAlignment = 'right';
Lbl_N_Elts.FontSize = 14;
Lbl_N_Elts.Layout.Row = 4;
Lbl_N_Elts.Layout.Column = 1;

 %Specify order of solution approximation: Langrange Shape/Interpolation Functions Order
  %Dropdown menu: lagrange shape function order
SF_order = uidropdown(g3);
SF_order.Items = {'Linear','Quadratic','Tertiary','Quartic','Quintic'};
SF_order.Layout.Row = 6;
SF_order.Layout.Column = 2;
SF_order.Value = 'Tertiary';
  %Label
Lbl_SF_order = uilabel(g3,'Text','Solution Order: ');
Lbl_SF_order.HorizontalAlignment = 'right';
Lbl_SF_order.FontSize = 14;
Lbl_SF_order.Layout.Row = 6;
Lbl_SF_order.Layout.Column = 1;

  %Dropdown menu: stabilization solving method
Method = uidropdown(g3);
Method.Items = {'Standard Galerkin Formulation','Isotropic Diffusion (Artificial Diffusion)'};
Method.Layout.Row = 8;
Method.Layout.Column = [2,3];
Method.Value = 'Standard Galerkin Formulation';
  %Label
Lbl_SF_order = uilabel(g3,'Text','Formulation: ');
Lbl_SF_order.HorizontalAlignment = 'right';
Lbl_SF_order.FontSize = 14;
Lbl_SF_order.Layout.Row = 8;
Lbl_SF_order.Layout.Column = 1;

  %Recalculate Button
recalc = uibutton(g3,'Text','Recalculate');
recalc.Layout.Row = 10;
recalc.Layout.Column = 4;
%=========================================

%Plotting and Heat Map Area
%
% C = C(Mesh); %Renumber with increasing distance between plates

% if Elts >= 34 %Reduce plotted points for better responsiveness
% post_pro_skip = floor(length(C)/(34*order));
% C = [C(1); C(2:post_pro_skip:end-1); C(end)];
% Nodes = length(C);
% end

% heat_map = heatmap(p_hm,C');

if Elts < 34
N = ceil(100/Elts); % of internal nodes evaluated and plotted, per Element
inter(order+1,N) = 0;
step = 2/(N+1);
for i=1:order+1
    inter(i,:) = polyval(psi(i,:),-1+step:step:1-step);
end

for i=1:Elts
y = [C(B(i,1)),C(B(i,:))'*inter,C(B(i,end))];
x = linspace((i-1)*h_e,i*h_e,N+2);
hold(C_plot,'on')
plot(C_plot,x,y,'b')
end

heat_map = heatmap(p_hm,C(Mesh)');

else % if Elts >= 34:  
   original_length = length(C);
   C = C(Mesh);
   x = linspace(0,1,original_length);
%    post_pro_skip = floor(original_length/(30*order));
   post_pro_skip = floor(original_length/(30*order+floor(Elts/1e3)));
   C_sample = [C(1); C(post_pro_skip:post_pro_skip:end-post_pro_skip); C(end)];
   x_sample = [x(1), x(post_pro_skip:post_pro_skip:end-post_pro_skip), x(end)];    

   plot(C_plot,x_sample,C_sample,'b')
   
%    C = [C(1); C(2:post_pro_skip:end-1); C(end)];
   Nodes = length(C_sample); 
   heat_map = heatmap(p_hm,C_sample');
end

heat_map.Title = "Heat Map of Cation Concentration";
heat_map.GridVisible = 'off';
heat_map.HandleVisibility = 'off';
heat_map.YDisplayLabels = '';
heat_map.Colormap = cool;
x_label([1,cast(Nodes/2,"int64"),Nodes],1) = ["0";num2str(L/2);num2str(L)];

% Adds more x-axis labels to the heat map
% total_x_labels = 4;
% for j=1:Nodes/(total_y_labels-1):Nodes
%    x_label(cast(j,"int64")) = num2str((j-1)/(Nodes-1));
% end

heat_map.XDisplayLabels = x_label;
heat_map.XLabel = "Distance {\itL} between plates";
%=========================================

%Callbacks
Width = 0.55*screen_size(1);
Height = 0.65*screen_size(2);
fui.CloseRequestFcn = @(src,event)closereq(src);
fui.Position = [(screen_size(1)-Width)/2 (screen_size(2)-Height)/2 Width Height];
fui.Resize = 'off';
% fui.AutoResizeChildren = 'off';
% fui.SizeChangedFcn = @(src,event)my_resize(src);
recalc.ButtonPushedFcn = @(src,event)recalculate(DD,uu,LL,N_Elts,SF_order,Method,C_plot,p_hm,g1);
%=========================================

%Set App Data
% setappdata(0,'D',DD);
% setappdata(0,'u',uu);
% setappdata(0,'L',LL);
% setappdata(0,'N_Elts',N_Elts);
% setappdata(0,'SF_order',SF_order);

toc
disp("GUI built. Displaying anytime!")
end






%FUNCTIONS
%=========
function closereq(f)
    selection = uiconfirm(f,'Are you sure you want to close this window?',...
            'Please Confirm');
        
    switch selection
        case 'OK'
            delete(f)
            diary off
            disp("GUI closed!")
        case 'Cancel'
            return
    end
end

% function my_resize(f)
%     dim = f.Position(3:4);
%     if dim(1) < 960
%         f.Resize = 'off';
%         f.Position(3) = 960;
%     end
%     if dim(2) < 600
%         f.Resize = 'off';
%         f.Position(4) = 600;
%     end
%     pause(0.1)
%     f.Resize = 'on';
%     f.AutoResizeChildren = 'on';
%     pause(0.1)
%     f.AutoResizeChildren = 'off';
% 
% end

function recalculate(DD,uu,LL,N_Elts,SF_order,Method,C_plot,p_hm,g1)
    param.u = uu.Value;
    param.L = LL.Value;
    param.Elts = N_Elts.Value;
    SF_order = SF_order.Value;
    
    selection = SF_order;
    switch selection
        case 'Linear'
            param.order = 1;
        case 'Quadratic'
            param.order = 2;
        case 'Tertiary'
            param.order = 3;
        case 'Quartic'
            param.order = 4;
        case 'Quintic'
            param.order = 5;
    end
    
%     delta = 0.5;
    selection = Method.Value;
    switch selection
        case 'Standard Galerkin Formulation'
            param.D = DD.Value;
        case 'Isotropic Diffusion (Artificial Diffusion)'
            param.D = DD.Value + param.u*(param.L/param.Elts)/(2*param.order);
    end
    
%Log
fprintf(1,"\n*\nRecalculating concentrations, with parameters:\n")
fprintf(1,"\t"+param.Elts+" Elements\n")
fprintf(1,"\tMaterial Diffusivity of "+param.D+" m^2/s\n")
fprintf(1,"\tFlow Velocity of "+param.u+" m/s\n")
fprintf(1,"\tDistance between the parallel plates of "+param.L+" meter(s)\n")
disp("Using:")
fprintf(1,"\t"+selection+"\n")
fprintf(1,"\tSolution approximated using Lagrange shape functions of order "+param.order+"\n")

%Get the Concentrations
tic
[C, B, psi, Mesh] = SolveForC(param);
disp("Recalculation Complete!")
toc
fprintf(1,"Updating GUI...\n")



%Update plot and heat map
tic
hold(C_plot,'off')
 
delete(p_hm);
p_hm = uipanel(g1,'Title','Heat Map');
p_hm.Layout.Row = 2;
p_hm.Layout.Column = 1;

h_e = param.L/param.Elts;
% N = cast(cast(150/param.Elts,"int64")+1,'like',1); % of internal nodes evaluated and plotted, per Element
if param.Elts < 34
    
N = ceil(100/param.Elts); % of internal nodes evaluated and plotted, per Element
inter(param.order+1,N) = 0;
step = 2/(N+1);

for i=1:param.order+1
    inter(i,:) = polyval(psi(i,:),-1+step:step:1-step);
end

for i=1:param.Elts
   y = [C(B(i,1)),C(B(i,:))'*inter,C(B(i,end))];
   x = linspace((i-1)*h_e,i*h_e,N+2);
   plot(C_plot,x,y,'b')
   hold(C_plot,'on') % To keep plots of previous elements and not delete them
end
 
heat_map = heatmap(p_hm,C(Mesh)');

else % for param.Elts >= 34:
    original_length = length(C);
%     post_pro_skip = floor(original_length/(30*param.order));
    post_pro_skip = floor(original_length/(30*param.order+floor(param.Elts/1e3)));
    C = C(Mesh);
    x_C = linspace(0,1,original_length);
    C = [C(1); C(post_pro_skip:post_pro_skip:end-post_pro_skip); C(end)];
    x_C = [x_C(1), x_C(post_pro_skip:post_pro_skip:end-post_pro_skip), x_C(end)];
    plot(C_plot,x_C,C,'b')
    
    heat_map = heatmap(p_hm,C');
end

%Common for any # of elements
heat_map.Title = "Heat Map of Cation Concentration";
heat_map.Colormap = cool;
heat_map.GridVisible = 'off';
heat_map.HandleVisibility = 'off';
heat_map.YDisplayLabels = '';
heat_map.XLabel = "Distance {\itL} between plates";
Nodes = length(C);
 
%clear labels
x_label(1:Nodes) = " ";
heat_map.XDisplayLabels = x_label;
%display new labels
x_label([1,cast(Nodes/2,"int64"),Nodes]) = ["0",num2str(param.L/2),num2str(param.L)];

heat_map.XDisplayLabels = x_label;

toc
disp("GUI updated. Displaying... Check GUI")
end