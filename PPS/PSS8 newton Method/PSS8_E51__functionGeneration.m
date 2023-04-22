clc;clear;

% Define the variables:
syms theta

% Define your function:
r = theta^3-1;

% Calculate the Jacobian:
dr = jacobian(r,theta);

% Save as a matlab function in a separate file:
matlabFunction(r,dr, 'File', 'fungen_r');
