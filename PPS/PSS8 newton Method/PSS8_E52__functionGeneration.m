clc;clear;

% Define the variables:
syms theta

% Define your function:
% phi = 0.5*(theta-1)^2;
phi = theta^2+5*exp(-theta^2)+theta;

% Calculate the Jacobian:
r  = jacobian(phi, theta);
dr = jacobian(r,theta) ;

% Save as a matlab function in a separate file:
matlabFunction(r,dr, 'file', 'fungen2_r','vars',{theta})
matlabFunction(phi,'file','fungen_phi_plot','vars',{theta})
