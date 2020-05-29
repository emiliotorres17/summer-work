% MAE471 - CFD
% Final Project
% Solving incompressible, 2D Navier-Stokes with Fractional Step Method on a
% Staggered Grid
 
%Method involves 3 steps:
%   I)      Solve Viscous Burgers Equation using FTCS
%   II)     Calculate Lagrange multiplier using Gauss-Seidel
%   III)    Project solution into subspace of solenoidal velocity fields
 
%%
 
clear all
clc
 
tic
 
%%%%%%%%%%%%    Inputs
Re  = 100; %Reynolds number
L   = 1; %length and height of cavity
nu  = 0.01; %kinematic viscosity
U   = Re*nu/L; %velocity at top wall
 
M   = 32; %number of cells in x
N   = 32; %number of cells in y
 
conv_crit = 10^-7; %convergence criteria
 
%%%%%%%%%%%%    Discretization
dx = L/M;
dy = L/N;
rh = 1/dx; %recipricol of h (dx = dy)
rh2 = 1/dx^2; %recipricol of h^2
rRe = 1/Re; %recipricol of Re
 
%%%%%%%%%%%%    Initial Conditions
u = zeros(M+1,N+2);
v = zeros(M+2,N+1);
u_star = u;
v_star = v;
phi = zeros(M+2,N+2);
res_gs = phi;
gs_RHS = phi;
gs_RHS2 = phi;
omega = zeros(M+1,N+1);
u_plot = zeros(1,N+1);
v_plot = zeros(1,N+1);
 
n = 0; %counter
t = 0; %time
n_gs_t = 0; %gauss-seidel counter (total number of iterations for all time)
conv = 1; %initial convergence value
 
%%%%%%%%%%%%    Boundary Conditions
%u velocity
u(:,N+2) = 2*U-u(:,N+1); %top wall
u(:,1) = -u(:,2); %bottom wall
u(1,:) = 0; %left wall
u(M+1,:) = 0; %right wall
%v velocity
v(:,N+1) = 0; %top wall
v(:,1) = 0; %bottom wall
v(1,:) = -v(2,:); %left wall
v(M+2,:) = -v(M+1,:); %right wall
 
%%%%%%%%%%%%    Calculate Solutions
 
while conv > conv_crit
    
    %%%%% find dt from stability checks
    %a & b values for constraints
    a_grid = u;
    b_grid = v;
    a = max(max(abs(a_grid)));
    b = max(max(abs(b_grid)));
 
    %following constraints, assume dx = dy
 
    %first stability constraint, parabolic part
    dt1 = 0.25*dx^2/nu;
 
    %second stability constraint, hyperbolic part 
    dt2 = dx/(a+b);
 
    %pick smallest dt
    if dt1 >= dt2
        dt = dt2;
    else
        dt = dt1;
    end
 
    %check third stability constraint, cell Reynolds number
    if (a+b)*dx/nu > 2
        disp('Unstable - Cell Reynolds Number')
    end
    
    %%%%% increase counters
    t = t+dt; %increase time by stable dt
    n = n+1; %counter
    t_plot(n) = t; %time for plotting
    n_gs = 0; %gauss-seidel counter for each set of loops
    
    %   Step I of Fractional Step Method
    
    %%%%% store solution at n for convergence check
    r_u = u;
    r_v = v;
    
    %%%%% find solutions to u @ time step n_star using FTCS
    for j = 2:N+1
        for i = 2:M
            %averages because of staggered mesh
            ua = 0.5*(u(i+1,j)+u(i,j)); %u_i+1,j
            ub = 0.5*(u(i,j)+u(i-1,j)); %u_i,j
            uc = 0.5*(u(i,j+1)+u(i,j)); %u_i+1/2,j+1/2
            ud = 0.5*(u(i,j-1)+u(i,j)); %u_i+1/2,j-1/2
    
            va = 0.5*(v(i+1,j)+v(i,j)); %v_i+1/2,j+1/2
            vb = 0.5*(v(i+1,j-1)+v(i,j-1)); %v_i+1/2,j-1/2
            
            %u-velocity star value
            u_star(i,j) = u(i,j)+dt*(-((ua^2-ub^2)*rh+(uc*va-ud*vb)*rh)...
                +rRe*((u(i+1,j)-2*u(i,j)+u(i-1,j))*rh2+(u(i,j+1)-2*u(i,j)...
                +u(i,j-1))*rh2));
        end
    end
    
    %%%%% find solutions to v @ time step n_star using FTCS
    for j = 2:N
        for i = 2:M+1
            %averages because of staggered mesh
            uc = 0.5*(u(i,j+1)+u(i,j)); %u_i+1/2,j+1/2
            ue = 0.5*(u(i-1,j)+u(i-1,j+1)); %u_i-1/2,j+1/2
            
            va = 0.5*(v(i+1,j)+v(i,j)); %v_i+1/2,j+1/2
            vc = 0.5*(v(i,j+1)+v(i,j)); %v_i,j+1
            vd = 0.5*(v(i,j)+v(i,j-1)); %v_i,j
            vf = 0.5*(v(i,j)+v(i-1,j)); %v_i-1/2,j+1/2
            
            %v-velocity star value
            v_star(i,j) = v(i,j)+dt*(-((uc*va-ue*vf)*rh+(vc^2-vd^2)*rh)...
                +rRe*((v(i+1,j)-2*v(i,j)+v(i-1,j))*rh2+(v(i,j+1)-2*v(i,j)+...
                v(i,j-1))*rh2));
        end
    end
    
    %%%%%%%%%%%%    Boundary Conditions
    %u_star velocity
    u_star(:,N+2) = 2*U-u_star(:,N+1); %top wall
    u_star(:,1) = -u_star(:,2); %bottom wall
    u_star(1,:) = 0; %left wall
    u_star(M+1,:) = 0; %right wall
    %v_star velocity
    v_star(:,N+1) = 0; %top wall
    v_star(:,1) = 0; %bottom wall
    v_star(1,:) = -v_star(2,:); %left wall
    v_star(M+2,:) = -v_star(M+1,:); %right wall
    
    %   Step II of Fractional Step Method
    
    %%%%% GS Convergence Condition
    %dynamically changing convergence criteria for Gauss-Seidel loop
    
    conv_gs = 1; %initial convergence value for Step II
    if conv*10^-2 >= 10^-3
        conv_gs_limit = 10^-3;
    elseif conv*10^-2 <= 10^-7
        conv_gs_limit = 10^-7;
    end
    
    %%%%% find solution to Langrange multiplier using Gauss-Seidel
    %iterates until Gauss-Seidel convergence criteria is met
    %pre-calculate right-hand sides of GS and residual equations for speed
    for j = 2:N+1
        for i = 2:M+1
            gs_RHS(i,j) = 0.25*(dx/dt)*(u_star(i,j)-u_star(i-1,j)+v_star(i,j)-v_star(i,j-1));
            gs_RHS2(i,j) = (1/dt)*(((u_star(i,j)-u_star(i-1,j))*rh)+((v_star(i,j)-v_star(i,j-1))*rh));
        end
    end
    
    while conv_gs > conv_gs_limit
        n_gs = n_gs+1; %gauss-seidel iteration counter
        n_gs_t = n_gs_t+1; %gauss-seidel iteration counter (total)
        
        %Gauss-Seidel
        for j = 2:N+1
            for i = 2:M+1
                phi(i,j) = 0.25*(phi(i-1,j)+phi(i+1,j)+phi(i,j-1)...
                    +phi(i,j+1))-gs_RHS(i,j);
            end
        end
        
        %%%%% update boundary conditions - Langrange multiplier ("pressure")
        phi(:,1) = phi(:,2); %bottom wall
        phi(:,N+2) = phi(:,N+1); %top wall
        phi(1,:) = phi(2,:); %left wall
        phi(M+2,:) = phi(M+1,:); %right wall  
 
    
        %%%%% check convergence of Lagrange multiplier
        for j = 2:N+1
            for i = 2:M+1
                %Residual for entire mesh
                res_gs(i,j) = (phi(i+1,j)-2*phi(i,j)+phi(i-1,j))*rh2 ...
                +(phi(i,j+1)-2*phi(i,j)+phi(i,j-1))*rh2 ...
                -gs_RHS2(i,j);
            end
        end
    
        conv_gs = max(max(abs(res_gs))); %check
    
    end
    
    %   Step III of Fractional Step Method
    
    %%%%% find u velocity solution at n+1
    for j = 2:N+1
        for i = 2:M
            u(i,j) = u_star(i,j)-(dt*rh)*(phi(i+1,j)-phi(i,j));
        end
    end
    
    %%%%% find v velocity solution at n+1
    for j = 2:N
        for i = 2:M+1
            v(i,j) = v_star(i,j)-(dt*rh)*(phi(i,j+1)-phi(i,j));
        end
    end
    
    %%%%% update boundary conditions - velocities
    %u velocity
    u(:,N+2) = 2*U-u(:,N+1); %top wall
    u(:,1) = -u(:,2); %bottom wall
    u(1,:) = 0; %left wall
    u(M+1,:) = 0; %right wall
    %v velocity
    v(:,N+1) = 0; %top wall
    v(:,1) = 0; %bottom wall
    v(1,:) = -v(2,:); %left wall
    v(M+2,:) = -v(M+1,:); %right wall
    
    %%%%% convergence of final solution
    res_u = (u-r_u)/dt; %du/dt
    res_v = (v-r_v)/dt; %dv/dt
    conv_u = max(max(abs(res_u)));
    conv_v = max(max(abs(res_v)));
    conv = max([conv_u conv_v]);
    
    %   Get Ready for Next Iteration
    
    %%%%% convergence history
    conv_hist_u(n) = conv_u;
    conv_hist_v(n) = conv_v;
 
    %%%%% Output Running Update
    fprintf('Convergence of u = %6.8f, Convergence of v = %6.8f\n',conv_u,conv_v);
    fprintf('Iteration %d, %d GS Iterations, %d Total GS Iterations\n',n,n_gs,n_gs_t);
end
 
toc
 
%%
 
%%%%%%%%%%%%    Plots
 
%%%%% vorticity
%calculate vorticity array
for j = 1:N+1
    for i = 1:M+1
        omega(i,j) = rh*(v(i+1,j)-v(i,j)-u(i,j+1)+u(i,j));
    end
end
 
%import ghia vorticity data
ghia_vor_6 = xlsread('ghia_vort.xlsx',1,'A3:B27');
ghia_vor_5 = xlsread('ghia_vort.xlsx',1,'C3:D25');
ghia_vor_4 = xlsread('ghia_vort.xlsx',1,'E3:F33');
ghia_vor_3 = xlsread('ghia_vort.xlsx',1,'G3:H39');
ghia_vor_2 = xlsread('ghia_vort.xlsx',1,'I3:J42');
ghia_vor_1 = xlsread('ghia_vort.xlsx',1,'K3:L41');
ghia_vor_0 = xlsread('ghia_vort.xlsx',1,'M3:N69');
ghia_vor_n1 = xlsread('ghia_vort.xlsx',1,'O3:P57');
ghia_vor_n2 = xlsread('ghia_vort.xlsx',1,'Q3:R37');
ghia_vor_n3 = xlsread('ghia_vort.xlsx',1,'S3:T31');
ghia_vor_n4 = xlsread('ghia_vort.xlsx',1,'U3:V28');
 
%plot vorticity
xy_vor = linspace(0,L,M+1);
con = [0.0 -0.5 0.5 -1.0 1.0 -2.0 2.0 -3.0 3.0 -4.0 -5.0]; %ghia vorticity values
figure
contour(xy_vor,xy_vor,omega',con) %plots contour at ghia values
hold on
%plot ghia values
plot(ghia_vor_6(:,1),ghia_vor_6(:,2),'ok')
plot(ghia_vor_5(:,1),ghia_vor_5(:,2),'ok')
plot(ghia_vor_4(:,1),ghia_vor_4(:,2),'ok')
plot(ghia_vor_3(:,1),ghia_vor_3(:,2),'ok')
plot(ghia_vor_2(:,1),ghia_vor_2(:,2),'ok')
plot(ghia_vor_1(:,1),ghia_vor_1(:,2),'ok')
plot(ghia_vor_0(:,1),ghia_vor_0(:,2),'ok')
plot(ghia_vor_n1(:,1),ghia_vor_n1(:,2),'ok')
plot(ghia_vor_n2(:,1),ghia_vor_n2(:,2),'ok')
plot(ghia_vor_n3(:,1),ghia_vor_n3(:,2),'ok')
plot(ghia_vor_n4(:,1),ghia_vor_n4(:,2),'ok')
%labeling
title(sprintf('Vorticity Contours for %d x %d Mesh',M,N),'fontsize',14)
xlabel('x','fontsize',14); ylabel('y','fontsize',14)
legend('Fractional Step Method','Ghia Values')
 
%%%%% velocities along lines through geometric center
x = linspace(0,L,M+2); %horizontal line
y = linspace(0,L,N+2); %vertical line
 
for q = 1:N+2
    u_plot(q) = 0.5*(u(N/2,q)+u((N+2)/2,q)); %u along vertical line
    v_plot(q) = 0.5*(v(q,(N+2)/2)+v(q,N/2)); %v along horizontal line
end
 
%ghia values for u and v along same lines
u_ghia = xlsread('ghia_u&v.xlsx',1,'A3:B19');
v_ghia = xlsread('ghia_u&v.xlsx',1,'D3:E19');
 
%plot velocities
figure
plot(y,u_plot,'-r'); hold on
plot(u_ghia(:,1),u_ghia(:,2),'ok');
title('u-Velocity Along Vertical Line Through Geometric Center of Cavity','fontsize',14);
xlabel('Vertical Line, y','fontsize',14); ylabel ('u-Velocity','fontsize',14)
legend(sprintf('%d x %d grid',M,N),'Ghia Values')
 
figure
plot(x,v_plot,'-r'); hold on
plot(v_ghia(:,1),v_ghia(:,2),'ok');
title('v-Velocity Along Horizontal Line Through Geometric Center of Cavity','fontsize',14);
xlabel('Horizontal Line, x','fontsize',14); ylabel('v-Velocity','fontsize',14)
legend(sprintf('%d x %d grid',M,N),'Ghia Values')
 
%plot convergence history for velocities
figure
plot(t_plot,log10(conv_hist_u),'-c'); hold on
plot(t_plot,log10(conv_hist_v),'-g');
title(sprintf('Convergence History for %d x %d Mesh',M,N),'fontsize',14)
xlabel('Time, s','fontsize',14); ylabel('Convergence, log scale','fontsize',14);
legend('u-velocity','v-velocity')
        1.2 GCI Accuracy Analysis
% MAE471 - CFD
% Final Project
% Solving incompressible, 2D Navier-Stokes for Lid-Driven Cavity
 
% Richard Extrapolation and GCI Analysis
% Determines accuracy for velocity solutions along geometric center lines
 
r = 2; %ratio between meshes, constant for 32 and 21
 
%u-velocity
f3u = 1.1089; %coarsest mesh, M = 32
f2u = 1.0523; %medium mesh, M = 64
f1u = 1.0258; %finest mesh, M = 128
 
p_u = log((f3u-f2u)/(f2u-f1u))/log(r); %order of convergence
f_u = f1u+((f1u-f2u)/(r^p_u-1)); %richardson extrapolation
 
GCI_12u = 1.25*(abs((f1u-f2u)/f1u)/(r^p_u-1));
GCI_23u = 1.25*(abs((f2u-f3u)/f2u)/(r^p_u-1));
 
fs_u = [f1u f2u f3u];
 
figure
plot([1 2 4],fs_u,'-or'); hold on
plot(0,f_u,'d');
title('Accuracy Analysis of u-Velocity Solution','fontsize',14);
ylabel('Maximum u-Velocity Along Vertical Center Line','fontsize',14);
xlabel('Normalized Grid Spacing','fontsize',14)
legend('Fractional Step Method','Richardson Extrapolation')
 
%v-velocity
f3v = 0.1675;
f2v = 0.1746;
f1v = 0.1775;
 
p_v = log((f3v-f2v)/(f2v-f1v))/log(r);
f_v = f1v+((f1v-f2v)/(r^p_v-1));
 
GCI_12v = 1.25*(abs((f1v-f2v)/f1v)/(r^p_v-1));
GCI_23v = 1.25*(abs((f2v-f3v)/f2v)/(r^p_v-1));
 
fs_v = [f1v f2v f3v];
 
figure
plot([1 2 4],fs_v,'-or'); hold on
plot(0,f_v,'d');
title('Accuracy Analysis of v-Velocity Solution','fontsize',14);
ylabel('Maximum v-Velocity Along Horizontal Center Line','fontsize',14);
xlabel('Normalized Grid Spacing','fontsize',14)
legend('Fractional Step Method','Richardson Extrapolation')
        1.3 ANSYS/Ghia Plotter
% MAE471 - CFD
% Final Project
% Solving incompressible, 2D Navier-Stokes in lid-drive cavity
 
%%Compares solution found using ANSYS Fluent to Ghia 1982
% 
%clear all
%clc
% 
%%%%%%   Import Data
%%import Ghia data
%u_ghia = xlsread('ghia_ANSYS.xlsx',1,'E4:F20');
%v_ghia = xlsread('ghia_ANSYS.xlsx',1,'B4:C20');
% 
%%import ANSYS data
%u_ANSYS = xlsread('ghia_ANSYS.xlsx',1,'E24:F66');
%v_ANSYS = xlsread('ghia_ANSYS.xlsx',1,'B24:C66');
% 
%%u-velocity plot
%figure
%plot(u_ANSYS(:,2),u_ANSYS(:,1),'-r'); hold on
%plot(u_ghia(:,1),u_ghia(:,2),'ok')
%title('u-Velocity Along Vertical Line Through Geometric Center of Cavity','fontsize',14);
%xlabel('Vertical Line, y','fontsize',14); ylabel ('u-Velocity','fontsize',14)
%legend('ANSYS Fluent','Ghia Values')
% 
%%v-velocity plot
%figure
%plot(v_ANSYS(:,1),v_ANSYS(:,2),'-r'); hold on
%plot(v_ghia(:,1),v_ghia(:,2),'ok')
%title('v-Velocity Along Horizontal Line Through Geometric Center of Cavity','fontsize',14);
%xlabel('Horizontal Line, x','fontsize',14); ylabel('v-Velocity','fontsize',14)
%legend('ANSYS Fluent','Ghia Values')
