%=========================================================================%
% Purpose:                                                                %
%   Final Project                                                         %
%   Solving incompressible, 2D Navier-Stokes with Fractional Step Method  %
%   on a Staggered Grid                                                   %
                                                                          %
% Method involves 3 steps:                                                %
%   I)      Solve Viscous Burgers Equation using FTCS                     %
%   II)     Calculate Lagrange multiplier using Gauss-Seidel              %
%   III)    Project solution into subspace of solenoidal velocity fields  %
%                                                                         %
% Author:                                                                 %
%   Emilio Torres                                                         %
%=========================================================================%
%-------------------------------------------------------------------------%
% Preamble                                                                %
%-------------------------------------------------------------------------%
close all
clear all
clc
tic
%-------------------------------------------------------------------------%
% Domain variables                                                        %
%-------------------------------------------------------------------------%
Re  = 100.0;                % Reynolds number
L   = 1;                    % length and height of cavity
nu  = 0.01;                 % kinematic viscosity
U   = Re*nu/L;              % velocity at top wall
M   = 32;                  % number of cells in x
N   = 32;                  % number of cells in y
%-------------------------------------------------------------------------%
% Steady state convergence criteria                                       %
%-------------------------------------------------------------------------%
conv_crit   = 10^-9;          % convergence criteria
tfinal      = 30.0;
%-------------------------------------------------------------------------%
% Discretization domain                                                   %
%-------------------------------------------------------------------------%
dx  = L/M;                              % x-step size
dy  = L/N;                              % y-step size
rh  = 1/dx;                             % reciprocal of h (dx = dy)
rh2 = 1/dx^2;                           % reciprocal of h^2
rRe = 1/Re;                             % reciprocal of Re
%-------------------------------------------------------------------------%
% Preallocating fields                                                    %
%-------------------------------------------------------------------------%
u       = zeros(M+1, N+2);              % u-velocity field
v       = zeros(M+2, N+1);              % v-velocity field
u_star  = zeros(M+2, N+2);              % u-star-velocity field (no pressure)             
v_star  = zeros(M+2, N+1);              % v-star-velocity field (no pressure)
phi     = zeros(M+2, N+2);              % Lagrangian multiplier
res_gs  = zeros(M+2, N+2);              % Guass Seidel residual          
gs_RHS  = zeros(M+2, N+2);              % right hand side of GS
gs_RHS2 = zeros(M+2, N+2);              % right hand side of GS
omega   = zeros(M+1, N+1);              % vorticity
u_plot  = zeros(1, N+1);                % center line u-velocity 
v_plot  = zeros(1, N+1);                % center line v-velocity
%-------------------------------------------------------------------------%
% Time counter variables                                                  %
%-------------------------------------------------------------------------%
n       = 0;        % counter
t       = 0.0;      % time
n_gs_t  = 0;        % gauss-seidel counter (total number of iterations for all time)
conv    = 1;        % initial convergence value
%-------------------------------------------------------------------------%
% Boundary conditions                                                     %
%-------------------------------------------------------------------------%
u(:,N+2)    = 2*U-u(:,N+1);         % u-velocity top wall
u(:,1)      = -u(:,2);              % u-velocity bottom wall
u(1,:)      = 0;                    % u-velocity left wall
u(M+1,:)    = 0;                    % u-velocity right wall
v(:,N+1)    = 0;                    % v-velocity top wall
v(:,1)      = 0;                    % v-velocity bottom wall
v(1,:)      = -v(2,:);              % v-velocity left wall
v(M+2,:)    = -v(M+1,:);            % v-velocity right wall
%=========================================================================%
% Time loop                                                               %
%=========================================================================%
while t < tfinal
    %---------------------------------------------------------------------%
    % time step calculation                                               %
    %---------------------------------------------------------------------%
    a_grid  = u;                            % a value for the grid
    b_grid  = v;                            % b value for the grid
    a       = max(max(abs(a_grid)));        % maximum a value
    b       = max(max(abs(b_grid)));        % maximum b value
    dt1     = 0.25*dx^2/nu;                 % parabolic time step constraint
    dt2     = dx/(a+b);                     % hyperbolic time step constraint
    dt      = min(dt1, dt2);                % calculating time step
    %---------------------------------------------------------------------%
    % Checking cell Reynolds number stability                             %
    %---------------------------------------------------------------------%
    if (a+b)*dx/nu > 2
        disp('Unstable - Cell Reynolds Number')
    end
    %---------------------------------------------------------------------%
    % Increasing counters and resetting counters                          %
    %---------------------------------------------------------------------%
    t           = t+dt;     % increase time by stable dt
    n           = n+1;      % counter
    t_plot(n)   = t;        % time for plotting
    n_gs        = 0;        % gauss-seidel counter for each set of loops
    %=====================================================================%
    % Step I of Fractional Step Method                                    %
    %=====================================================================%
    %---------------------------------------------------------------------%
    % Storing solutions to time step t^{n}                                %
    %---------------------------------------------------------------------%
    r_u = u;                % u-vel. solution for t^{n} for convergence
    r_v = v;                % v-vel. solution for t^{n} for convergence
    %---------------------------------------------------------------------%
    % Calculating u-star values                                           %
    %---------------------------------------------------------------------%
    for j = 2:N+1
        for i = 2:M
            %-------------------------------------------------------------%
            % Calculating u-vel. and v-vel. average values                %
            %-------------------------------------------------------------%
            ua = 0.5*(u(i+1,j)+u(i,j));             % u_i+1,j
            ub = 0.5*(u(i,j)+u(i-1,j));             % u_i,j
            uc = 0.5*(u(i,j+1)+u(i,j));             % u_i+1/2,j+1/2
            ud = 0.5*(u(i,j-1)+u(i,j));             % u_i+1/2,j-1/2
            va = 0.5*(v(i+1,j)+v(i,j));             % v_i+1/2,j+1/2
            vb = 0.5*(v(i+1,j-1)+v(i,j-1));         % v_i+1/2,j-1/2
            %-------------------------------------------------------------%
            % u-star values                                               %
            %-------------------------------------------------------------%
            u_star(i,j) = u(i,j)+dt*(-((ua^2-ub^2)*rh+(uc*va-ud*vb)*rh)...
                +rRe*((u(i+1,j)-2*u(i,j)+u(i-1,j))*rh2+(u(i,j+1)-2*u(i,j)...
                +u(i,j-1))*rh2));
        end
    end
    %---------------------------------------------------------------------%
    % Calculating v-star values                                           %
    %---------------------------------------------------------------------%
    for j = 2:N
        for i = 2:M+1
            %-------------------------------------------------------------%
            % Calculating u-vel. and v-vel. average values                %
            %-------------------------------------------------------------%
            uc = 0.5*(u(i,j+1)+u(i,j));         %u_i+1/2,j+1/2
            ue = 0.5*(u(i-1,j)+u(i-1,j+1));     %u_i-1/2,j+1/2
            va = 0.5*(v(i+1,j)+v(i,j));         %v_i+1/2,j+1/2
            vc = 0.5*(v(i,j+1)+v(i,j));         %v_i,j+1
            vd = 0.5*(v(i,j)+v(i,j-1));         %v_i,j
            vf = 0.5*(v(i,j)+v(i-1,j));         %v_i-1/2,j+1/2
            %-------------------------------------------------------------%
            % Calculating v-star values                                   %
            %-------------------------------------------------------------%
            v_star(i,j) = v(i,j)+dt*(-((uc*va-ue*vf)*rh+(vc^2-vd^2)*rh)...
                +rRe*((v(i+1,j)-2*v(i,j)+v(i-1,j))*rh2+(v(i,j+1)-2*v(i,j)+...
                v(i,j-1))*rh2));
        end
    end
    %---------------------------------------------------------------------%
    % Boundary conditions                                                 %
    %---------------------------------------------------------------------%
    u_star(:,N+2)   = 2*U-u_star(:,N+1);    % u-star top wall
    u_star(:,1)     = -u_star(:,2);         % u-star bottom wall
    u_star(1,:)     = 0;                    % u-star left wall
    u_star(M+1,:)   = 0;                    % u-star right wall
    v_star(:,N+1)   = 0;                    % v-star top wall
    v_star(:,1)     = 0;                    % v-star bottom wall
    v_star(1,:)     = -v_star(2,:);         % v-star left wall
    v_star(M+2,:)   = -v_star(M+1,:);       % v-star right wall
    %=====================================================================%
    % Step II of Fractional Step Method                                   %
    %=====================================================================%
    %---------------------------------------------------------------------%
    % Dynamic GS convergence conditions                                   %
    %---------------------------------------------------------------------%
    conv_gs = 1;                    % resetting the GS convergence
    if conv*10^-2 >= 10^-3 
        conv_gs_limit = 10^-9;      % convergence criteria for early in time
    elseif conv*10^-2 <= 10^-7
        conv_gs_limit = 10^-9;      % convergence criteria near steady state
    end
    %---------------------------------------------------------------------%
    % Finding solution to Lagrange multiplier using Gauss-Seidel          %
    % iterates until Gauss-Seidel convergence criteria is met             %
    %---------------------------------------------------------------------%
    for j = 2:N+1       % pre-calculate RHS of GS and residual equations 
        for i = 2:M+1
            gs_RHS(i,j) = 0.25*(dx/dt)*(u_star(i,j)-u_star(i-1,j)+ ...
                            v_star(i,j)-v_star(i,j-1));
            gs_RHS2(i,j) = (1/dt)*(((u_star(i,j)-u_star(i-1,j))*rh)+ ...
                            ((v_star(i,j)-v_star(i,j-1))*rh));
        end
    end                 
    %---------------------------------------------------------------------%
    % GS iterations                                                       %
    %---------------------------------------------------------------------%
    while conv_gs > conv_gs_limit
        n_gs = n_gs+1;                      % GS iteration counter
        n_gs_t = n_gs_t+1;                  % GS iteration counter (total)
        for j = 2:N+1
            for i = 2:M+1
                phi(i,j) = 0.25*(phi(i-1,j)+phi(i+1,j)+phi(i,j-1)...
                    +phi(i,j+1))-gs_RHS(i,j);
            end
        end
        %-----------------------------------------------------------------%
        % Update Lagrange multiplier ("pressure") BCs                     %
        %-----------------------------------------------------------------%
        phi(:,1)    = phi(:,2);             % phi bottom wall
        phi(:,N+2)  = phi(:,N+1);           % phi top wall
        phi(1,:)    = phi(2,:);             % phi left wall
        phi(M+2,:)  = phi(M+1,:);           % phi right wall  
        %-----------------------------------------------------------------%
        % Check convergence of Lagrange multiplier                        %
        %-----------------------------------------------------------------%
        for j = 2:N+1
            for i = 2:M+1
                res_gs(i,j) = (phi(i+1,j)-2*phi(i,j)+phi(i-1,j))*rh2 ...
                                +(phi(i,j+1)-2*phi(i,j)+phi(i,j-1))*rh2 ...
                                -gs_RHS2(i,j);  % residual for entire mesh
            end
        end
        conv_gs = max(max(abs(res_gs)));        % check convergence
    end
    %=====================================================================%
    % Step III of Fractional Step Method                                  %
    %=====================================================================%
    %---------------------------------------------------------------------%
    % Calculating u-velocity solution at n+1                              %
    %---------------------------------------------------------------------%
    for j = 2:N+1
        for i = 2:M
            u(i,j) = u_star(i,j)-(dt*rh)*(phi(i+1,j)-phi(i,j));
        end
    end
    %---------------------------------------------------------------------%
    % Calculating v-velocity solution at n+1                              %
    %---------------------------------------------------------------------%
    for j = 2:N
        for i = 2:M+1
            v(i,j) = v_star(i,j)-(dt*rh)*(phi(i,j+1)-phi(i,j));
        end
    end
    %---------------------------------------------------------------------%
    % Update velocities BCs                                               %
    %---------------------------------------------------------------------%
    u(:,N+2)    = 2*U-u(:,N+1);     % u-velocity top wall
    u(:,1)      = -u(:,2);          % u-velocity bottom wall
    u(1,:)      = 0;                % u-velocity left wall
    u(M+1,:)    = 0;                % u-velocity right wall
    v(:,N+1)    = 0;                % v-velocity top wall
    v(:,1)      = 0;                % v-velocity bottom wall
    v(1,:)      = -v(2,:);          % v-velocity left wall
    v(M+2,:)    = -v(M+1,:);        % v-velocity right wall
    %---------------------------------------------------------------------%
    % Checking the velocity convergence                                   %
    %---------------------------------------------------------------------%
    res_u   = (u-r_u)/dt;                   % rate change of u-velocity
    res_v   = (v-r_v)/dt;                   % rate change of v-velocity
    conv_u  = max(max(abs(res_u)));         % u-vel. convergence
    conv_v  = max(max(abs(res_v)));         % v-vel. convergence
    conv    = max([conv_u conv_v]);         % maximum vel. rate of change
    %---------------------------------------------------------------------%
    % Storing convergence history                                         %
    %---------------------------------------------------------------------%
    conv_hist_u(n) = conv_u;                % u-velocity convergence
    conv_hist_v(n) = conv_v;                % v-velocity convergence
    %---------------------------------------------------------------------%
    % Time step print statement                                           %
    %---------------------------------------------------------------------%
    fprintf('time = %10.5e\n', t);
    fprintf('Convergence of u = %6.8e, Convergence of v = %6.8e\n' ...
                ,conv_u,conv_v);
    fprintf('Iteration %d, %d GS Iterations, %d Total GS Iterations\n' ...
                ,n,n_gs,n_gs_t);
end
toc
%=========================================================================%
% Post processing                                                         %
%=========================================================================%
%-------------------------------------------------------------------------%
% Calculating vorticity                                                   %
%-------------------------------------------------------------------------%
for j = 1:N+1
    for i = 1:M+1
        omega(i,j) = rh*(v(i+1,j)-v(i,j)-u(i,j+1)+u(i,j));
    end
end
%-------------------------------------------------------------------------%
% Plotting vorticity                                                      %
%-------------------------------------------------------------------------%
xy_vor      = linspace(0, L, M+1);              % XY for plotting
con         = [0.0 -0.5 0.5 -1.0 1.0 -2.0 2.0 -3.0 ...
                    3.0 -4.0 -5.0];             % vorticity color levels
figure
contour(xy_vor,xy_vor,omega',con)               % vorticity
title(sprintf('Vorticity Contours for %d x %d Mesh',M,N),'fontsize',14,...
            'interpreter', 'latex');
xlabel('x','fontsize',14, 'interpreter', 'latex'); 
ylabel('y','fontsize',14, 'interpreter', 'latex');
legend('Fractional Step Method', 'interpreter', 'latex');
%-------------------------------------------------------------------------%
% Plotting velocity on the geometric center                               %
%-------------------------------------------------------------------------%
x = linspace(0,L,M+2);                          % horizontal center line
y = linspace(0,L,N+2);                          % vertical center line
for q = 1:N+2
    u_plot(q) = 0.5*(u(N/2,q)+u((N+2)/2,q));    % u along vertical line
    v_plot(q) = 0.5*(v(q,(N+2)/2)+v(q,N/2));    % v along horizontal line
end
%-------------------------------------------------------------------------%
% Plotting u-velocity on the geometric center                             %
%-------------------------------------------------------------------------%
figure
hold on
plot(y,u_plot,'-r');
title('u-Velocity Along Vertical Line Through Geometric Center of Cavity',...
            'fontsize',14, 'interpreter', 'latex');
xlabel('Vertical Line, y','fontsize',14, 'interpreter', 'latex'); 
ylabel ('u-Velocity','fontsize',14, 'interpreter', 'latex');
legend(sprintf('%d x %d grid',M,N), 'interpreter', 'latex');
%-------------------------------------------------------------------------%
% Plotting v-velocity on the geometric center                             %
%-------------------------------------------------------------------------%
figure
hold on
plot(x,v_plot,'-r'); 
title('v-Velocity Along Horizontal Line Through Geometric Center of Cavity', ...
            'fontsize',14, 'interpreter', 'latex');
xlabel('Horizontal Line, x','fontsize',14, 'interpreter', 'latex'); 
ylabel('v-Velocity','fontsize',14, 'interpreter', 'latex');
legend(sprintf('%d x %d grid',M,N), 'interpreter', 'latex');
%-------------------------------------------------------------------------%
% Plotting v-velocity on the geometric center                             %
%-------------------------------------------------------------------------%
figure
hold on
plot(t_plot,log10(conv_hist_u),'-c'); 
plot(t_plot,log10(conv_hist_v),'-g');
title(sprintf('Convergence History for %d x %d Mesh',M,N),'fontsize',14, ...
            'interpreter', 'latex');
xlabel('Time, s','fontsize',14, 'interpreter', 'latex'); 
ylabel('Convergence, log scale','fontsize',14, 'interpreter', 'latex');
legend('u-velocity','v-velocity', 'interpreter', 'latex');
%=========================================================================%
% GCI Accuracy Analysis                                                   %
%=========================================================================%
%-------------------------------------------------------------------------%
% u-velocity analysis                                                     %
%-------------------------------------------------------------------------%
r       = 2;                                        % ratio between meshes
f3u     = 1.1089;                                   % u-velocity coarsest mesh, M = 32
f2u     = 1.0523;                                   % u-velocity medium mesh, M = 64
f1u     = 1.0258;                                   % u-velocity finest mesh, M = 128
p_u     = log((f3u-f2u)/(f2u-f1u))/log(r);          % v-velocity order of convergence
f_u     = f1u+((f1u-f2u)/(r^p_u-1));                % Richardson extrapolation
GCI_12u = 1.25*(abs((f1u-f2u)/f1u)/(r^p_u-1));      % GCI value 1-2 
GCI_23u = 1.25*(abs((f2u-f3u)/f2u)/(r^p_u-1));      % GCI value 2-3
fs_u    = [f1u f2u f3u];                            % u-vel. solution array
%-------------------------------------------------------------------------%
% u-velocity GCI plots                                                    %
%-------------------------------------------------------------------------%
figure
hold on
plot([1 2 4],fs_u,'-or'); 
plot(0,f_u,'d');
title('Accuracy Analysis of u-Velocity Solution','fontsize',14, 'interpreter', 'latex');
ylabel('Maximum u-Velocity Along Vertical Center Line','fontsize',14, 'interpreter', 'latex');
xlabel('Normalized Grid Spacing','fontsize',14, 'interpreter', 'latex');
legend('Fractional Step Method','Richardson Extrapolation', 'interpreter', 'latex');
%-------------------------------------------------------------------------%
% v-velocity GCI plots                                                    %
%-------------------------------------------------------------------------%
f3v     = 0.1675;                                   % v-velocity coarsest mesh, M = 32
f2v     = 0.1746;                                   % v-velocity medium mesh, M = 64
f1v     = 0.1775;                                   % velocity finest mesh, M = 128
p_v     = log((f3v-f2v)/(f2v-f1v))/log(r);          % v-velocity order of convergence 
f_v     = f1v+((f1v-f2v)/(r^p_v-1));                % v-velocity Richardson extrapolation
GCI_12v = 1.25*(abs((f1v-f2v)/f1v)/(r^p_v-1));      % GCI value 1-2
GCI_23v = 1.25*(abs((f2v-f3v)/f2v)/(r^p_v-1));      % GCI value 2-3
fs_v    = [f1v f2v f3v];                            % v-vel. solution array 
%-------------------------------------------------------------------------%
% v-velocity GCI plots                                                    %
%-------------------------------------------------------------------------%
figure
plot([1 2 4],fs_v,'-or'); hold on
plot(0,f_v,'d');
title('Accuracy Analysis of v-Velocity Solution','fontsize',14, 'interpreter', 'latex');
ylabel('Maximum v-Velocity Along Horizontal Center Line','fontsize',14, 'interpreter', 'latex');
xlabel('Normalized Grid Spacing','fontsize',14, 'interpreter', 'latex')
legend('Fractional Step Method','Richardson Extrapolation', 'interpreter', 'latex')
