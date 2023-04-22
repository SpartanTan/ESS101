clc;clear;
close all

%% Defining the function: (It is better to define this in a separate file)
n = 2;
x  = sym('x',[n,1]);
syms t real

f = [ x(2);
      5*(1-x(1)^2)*x(2) - x(1)];

matlabFunction(f, 'file', 'VanDerPol','vars',{t,x});
clear x t f

%% Defining r(xNext,x,u): (It is better to define this in a separate file)
xNext = sym('xNext',[2,1]);
x = sym('x',[2,1]);
syms t dt real

r = x + dt * VanDerPol(t+dt,xNext) - xNext;
dr = jacobian(r,xNext);

matlabFunction(r,dr, 'file', 'rFileIEuler','vars',{t,x,xNext,dt});
clear xNext x t dt r dr

%% Simulation parameters:
tFinal = 10;
dtEuler = 0.1;
x0 = [1;0];

%% The Solution:
% solve using ode15s and store in t.ode15s and x.ode15s
[t.ode15s,x.ode15s] = ode15s(@VanDerPol, [0 tFinal], x0);

% Plot the results:
% figure(1);clf
% for subfig = 1:n
%     subplot(n,1,subfig);hold on
%     plot(t.ode15s,x.ode15s(:,subfig),'k','marker','.','markersize',10)
% end
% 
% disp(['Please press any key to continue...'])
% pause

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
        
        xNext = xNext - dr\r;
        
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
% figure(2);clf
clf
for subfig = 1:n
    subplot(n,1,subfig);hold on
    plot(t.ode15s,x.ode15s(:,subfig),'k','marker','.','markersize',10)
    plot(t.euler,x.euler(subfig,:),'r')
    legend('ode15s','Implicit Euler');
end
