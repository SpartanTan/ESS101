% PSS 9 exercise 7.4 (Euler)
clc;clear;
close all

%% Defining the function:
f = @(t,x)(x/5+10*exp(t/5)*cos(10*t));

%% Simulation parameters:
tFinal = 3;
dtEuler = 0.001;
dtRK2 = 0.01;

%% The exact Solution:
tExact = 0:0.01:tFinal;
xExact = exp(tExact/5).*sin(10*tExact);

%% Explicit Euler Method:
nEuler = tFinal / dtEuler;

t = zeros(nEuler,1);
xEuler = zeros(nEuler,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Please write your code here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:length(t)
    xEuler(i) = xEuler(i-1) + dtEuler * f(t(i-1),xEuler(i-1));
    t(i) = t(i-1)+0.001;
end
%% Runge Kutta Method of order 2:
a=1;
b1 = 1/2;
b2 = 1/2;

nRK2 = tFinal / dtRK2;

t = zeros(nRK2,1);
xRK2 = zeros(nRK2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Please write your code here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:nRK2
    K1 = f(t(i-1), xRK2(i-1));
    K2 = f(t(i-1), xRK2(i-1)+dtRK2*K1);
    xRK2(i) = xRK2(i-1) + dtRK2*(b1*K2+b2*K1);
    t(i) = t(i-1) + dtRK2;
end
%% Compare the results:
plot(tExact,xExact  ,'marker','.','markersize',20)
hold on
plot(t,xEuler  ,'marker','.','markersize',10)
plot(t,xRK2      ,'marker','.','markersize',10)
xlabel('t');
ylabel('x');
legend('Exact','Euler Method','RK2')

%% Global error:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Please write your code here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_Euler = zeros(N,1);
err_RK2 = zeros(N,1); 
N_table = floor(logspace(1,3,50));
for n=1:length(N_table)
    N = N_table(n);
    dt = tFinal/N_table(n);
    t = zeros(N,1);
    xEuler = zeros(N,1);
    a=1;
    b1 = 1/2;
    b2 = 1/2;
    
    nRK2 = tFinal / dtRK2;
    xRK2 = zeros(N,1);
    for i=2:length(t)
        xEuler(i) = xEuler(i-1) + dt * f(t(i-1),xEuler(i-1));
        t(i) = t(i-1)+dt;

        K1 = f(t(i-1), xRK2(i-1));
        K2 = f(t(i-1), xRK2(i-1)+dtRK2*K1);
        xRK2(i) = xRK2(i-1) + dt*(b1*K2+b2*K1);
        t(i) = t(i-1) + dt;
    end
    err_Euler(n) = norm(xEuler(end)-xExact(end),inf);
    err_RK2(n) = norm(xRK2(end) - xExact(end), "inf");
end

plot(dt_list, err_RK2)
hold on
plot(dt_list, err_Euler)