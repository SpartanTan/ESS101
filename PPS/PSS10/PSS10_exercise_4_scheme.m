% PSS 9 exercise 7.4 (Euler)
clc;clear;
close all

%% Defining the function:
f = @(t,x)(x/5+10*exp(t/5)*cos(10*t));

%% Simulation parameters:
tFinal  = 3;
dtEuler = 0.2;
dtRK2   = 0.2;

%% The exact Solution:
tExact = 0:0.01:tFinal;
xExact = exp(tExact/5).*sin(10*tExact);

%% Explicit Euler Method:

nEuler = tFinal / dtEuler;

tEuler = zeros(nEuler,1);
xEuler = zeros(nEuler,1);

for j = 2:nEuler+1
    tEuler(j) = tEuler(j-1) + dtEuler;
    xEuler(j) = xEuler(j-1) + dtEuler*f(tEuler(j-1),xEuler(j-1));
end

%% Runge Kutta Method of order 2:
a=1;
b1 = 1/2;
b2 = 1/2;

nRK2 = tFinal / dtRK2;

tRK2 = zeros(nRK2,1);
xRK2 = zeros(nRK2,1);

for j = 2:nRK2+1
    tRK2(j) = tRK2(j-1) + dtRK2;
    K1 = f( tRK2(j-1) , xRK2(j-1) );
    K2 = f( tRK2(j-1)+a*dtRK2 , xRK2(j-1)+a*dtRK2*K1);
    xRK2(j) = xRK2(j-1) + dtRK2 * b1 * K1 + dtRK2 * b2 * K2 ;
end

%% Compare the results:
plot(tExact,xExact  ,'marker','.','markersize',20)
hold on
plot(tEuler,xEuler  ,'marker','.','markersize',10)
plot(tRK2,xRK2      ,'marker','.','markersize',10)
xlabel('t');
ylabel('x');
legend('Exact','Euler Method','RK2')
