%=========================================================================%
% Purpose:                                                                %
%   The purpose of this script is to generate the vorticity contour plots %
%   for the various meshes.                                               %
%                                                                         %
% Author:                                                                 %
%   Emilio Torres                                                         %
%=========================================================================%                                                                      
close all
clear all
clc
%-------------------------------------------------------------------------%
% Plotting vorticity                                                      %
%-------------------------------------------------------------------------%
M   = [32, 64, 128, 256];
for i = 1:length(M)
    fname       = strcat(['../../data/2d-solutions/vorticity-', num2str(M(i)), '.dat']);
    omega       = importdata(fname);
    xy_vor      = linspace(0, 1.0, M(i)+1);              % XY for plotting
    con         = [0.0 -0.5 0.5 -1.0 1.0 -2.0 2.0 -3.0 ...
                        3.0 -4.0 -5.0];             % vorticity color levels
    figure
    contour(xy_vor,xy_vor,omega',con)               % vorticity
    title(sprintf('Vorticity Contours for %d x %d Mesh',M(i),M(i)),'fontsize',14,...
                'interpreter', 'latex');
    xlabel('x','fontsize',14, 'interpreter', 'latex'); 
    ylabel('y','fontsize',14, 'interpreter', 'latex');
    legend('Fractional Step Method', 'interpreter', 'latex');
    name    = strcat(['vorticity-', num2str(M(i))]);
    saveas(gcf, strcat(['../media/', name, '.png']))
end
