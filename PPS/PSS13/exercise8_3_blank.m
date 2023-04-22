% This script contains an implicit Euler integrator and applies it to the 
% Van der Pol oscillator.
% There are 4 lines for the students to fill.

clc;clear;
close all

%% Defining the function: (It is better to define this in a separate file)
n = 2;
x  = sym('x',[n,1]);
syms t real

%%%%%%%%%%%%%%%%%%%%%%
% f = ...; % FILL HERE
%%%%%%%%%%%%%%%%%%%%%%

matlabFunction(f, 'file', 'VanDerPol','vars',{t,x});
clear x t f

%% Defining r(xNext,x,u): (It is better to define this in a separate file)
xNext = sym('xNext',[2,1]);
x = sym('x',[2,1]);
syms t dt real

%%%%%%%%%%%%%%%%%%%%%%%
% r = ...;  % FILL HERE
% dr = ...; % FILL HERE
%%%%%%%%%%%%%%%%%%%%%%%

matlabFunction(r,dr, 'file', 'rFileIEuler','vars',{t,x,xNext,dt});
clear xNext x t dt r dr

%% Simulation parameters:
tFinal = 3;
dtEuler = 0.1;
x0 = [1;0];

%% The Solution:
% solve using ode15s and store in t.ode15s and x.ode15s
[t.ode15s,x.ode15s] = ode15s(@VanDerPol, [0 tFinal], x0);

% Plot the results:
figure(1);clf
for subfig = 1:n
    subplot(n,1,subfig);hold on
    plot(t.ode15s,x.ode15s(:,subfig),'k','marker','.','markersize',10)
end

disp(['Please press any key to continue...'])
pause

%% Simulation:
% Simulate using IEuler
Nsteps = tFinal/dtEuler;
t.euler = [0:dtEuler:tFinal];
x.euler = [x0,zeros(n,Nsteps-1)];
xNext = x0;

% Loop for the Implicit Euler
for k = 1:Nsteps
    % Newton iteration
    iter = true;
    alpha = 1;
    niter = 0;
    while iter
        [r,dr] = rFileIEuler(t,x.euler(:,k),xNext,dtEuler);
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % xNext = ; % FILL HERE
        %%%%%%%%%%%%%%%%%%%%%%%
        
        norm(r);
        if norm(r) < 1e-5
            iter = false;
        else
            niter = niter + 1;
        end
    end
    x.euler(:,k+1) = xNext;
end

% Plot the results:
figure(2);clf
for subfig = 1:n
    subplot(n,1,subfig);hold on
    plot(t.euler,x.euler(subfig,:),'r')
end


