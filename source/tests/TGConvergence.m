%% Plot the convergence of the TG Vortex

close all
clear all
clc

% Spatial
data = importdata('TGVortexErrorOneStep.txt');

N = data(11:end,1).^0.5;
error = data(11:end,2)*100;

figure; hold on;
loglog(N, error, 'o-', 'LineWidth', 3)
xlabel('$\sqrt{N_{CVs}}$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$Error \ [\%]$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 15)
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
box on;
xx = [50 500];
yy = [50 5];
loglog(xx, yy, '-', 'LineWidth', 2)
L = legend('$Error$', '$1^{st} \ Order$');
set(L, 'interpreter', 'latex', 'fontsize', 15)


% Temporal
data = importdata('TGVortexErrorTime.txt');

dt = data(1:end,1);
error = data(1:end,2)*100;

figure; hold on;
plot(dt, error, 'o-', 'LineWidth', 3)
xlabel('$\Delta t \ [s]$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$Error \ [\%]$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 15)
set(gca, 'xscale', 'log')
%set(gca, 'yscale', 'log')
box on;
% xx = [50 500];
% yy = [50 5];
% loglog(xx, yy, '-', 'LineWidth', 2)
L = legend('$Error$', '$1^{st} \ Order$');
set(L, 'interpreter', 'latex', 'fontsize', 15)

