%% AE 370 Abaqus Project
% Mike O'Connell, Max Feinberg
% 
% Convergence analysis for part (c)
function ae370abaqus
clc;close all;
%% Import Data
%       [element size(ul), max end point deflection(mm), number of nodes]
% data from 3-node plane stress triangular elements
data_a =   [20, .1687, 14;
            10, .5449, 24;
            5, 1.993, 46;
            3.75, 3.207, 60;
            2.5, 6.030, 90;
            1.75, 9.226, 128;
            1, 13.97, 222;
            0.5, 33.48, 663;
            0.25, 51.48, 2205;
            .125, 59.48, 7929;
            0, 60.15, 10010;
            0, 60.61, 12111;
            0, 60.98, 14676;
            0, 61.43, 19264];
        
% data from 6-node plane stress triangular elements
data_b =   [20, 39.23, 39;
            10, 46.13, 69;
            5, 51.143, 135;
            3.75, 52.70, 177;
            2.5, 53.97, 267;
            1.75, 54.80, 381;
            1, 56.20, 663;
            0.5, 60.96, 2205;
            0.25, 62.28, 7929;
            0.155, 62.54, 18473];
        
%% Test for convergence
% I will test for convergence by showing that the condition number
%
% First, get the delta u_i = u_i - u_i-1
conditionnum_a = convergence_test(data_a)
conditionnum_b = convergence_test(data_b)

% Test for convergence by checking when the condition number falls below 
% .05. Thus, furhter increasing the mesh density will not significantly 
% change the end point deflection. This suggests that we have converged 
% close to the exact solution.  Looking at the outputs from the above
% lines, we see that both a and b converge on the solution.

% Note: the convergence number was used instead of the logarithmic
% derivative because it accounts for the change in step size between data
% points.


%% Make pretty pictures
figure(1)
hold on;
plotter(data_a,'3-node plane stress','bx-');
title_string = {'\fontsize{12}Maximum End Point Deflection vs Mesh Density'
    '\fontsize{12}Using Triangular Elements and 3-Node Plane Stress'
    '\fontsize{7}AE 370 Abaqus Finite Element Project'
    '\fontsize{7}Max Feinberg, Mike O''Connell'};
title(title_string)

figure(2)
hold on;
plotter(data_b,'6-node plane stress','ro-');
title_string(2) = {'\fontsize{12}Using Triangular Elements and 6-Node Plane Stress'};
title(title_string)

figure(3)
hold on;
plotter(data_a,'3-node plane stress','bx-');
plotter(data_b,'6-node plane stress','ro-');
legend('show','Location','east')
title_string(2) = {'\fontsize{12}Using Triangular Elements'};
title(title_string)

% Save it all
fprintf('%s','saving figure (1) to ae370abaqus_3node.png...')
saveas(figure(1),'ae370abaqus_3node.png')
fprintf('%s\n%s','done','saving figure (2) to ae370abaqus_6node.png...')
saveas(figure(2),'ae370abaqus_6node.png')
fprintf('%s\n%s','done','saving figure (3) to ae370abaqus_combined.png...')
saveas(figure(3),'ae370abaqus_combined.png')
fprintf('%s\n\n','done')

function [handle] = plotter(data,displayname,line_options)
handle = plot(data(:,3)./110,data(:,2),line_options,'DisplayName',displayname);
xlabel('mesh density (nodes/m^2)')
ylabel('end point deflection (mm)')
set(gca,'xscale','log')

% Test convergence by estimating condition number at point of interest.
function [retval] = convergence_test(data)
retval = data(2:end,3)./(data(2:end,3)-data(1:end-1,3)).*(data(2:end,2)-data(1:end-1,2))./data(2:end,2);





















