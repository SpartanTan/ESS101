% PSS 9 exercise 7.4 (Euler)

clc;clear;
close all

%% Defining the function:
f = @(t,x)(x/5+10*exp(t/5)*cos(10*t));

%% Simulation parameters:
tFinal  = 3;
a       = 1;
b1      = 1/2;
b2      = 1/2;

%% The exact Solution:
tExact = 0:0.01:tFinal;
xExact = exp(tExact/5).*sin(10*tExact);


%% Integration for various time steps:
N_table = floor(logspace(1,3,50));

for index = 1:length(N_table)
    
    N = N_table(index);
    dt(index) = tFinal / N_table(index);
    t = zeros(N,1);
    xEuler = zeros(N,1);
    xRK2 = zeros(N,1);
    
    for j = 2:N+1
        
        t(j) = t(j-1) + dt(index);
        
        % Euler step
        xEuler(j) = xEuler(j-1) + dt(index)*f(t(j-1),xEuler(j-1));
        
        % RK2 step
        K1 = f( t(j-1) , xRK2(j-1) );
        K2 = f( t(j-1)+a*dt(index) , xRK2(j-1)+a*dt(index)*K1);
        xRK2(j) = xRK2(j-1) + dt(index) * b1 * K1 + dt(index) * b2 * K2 ;

    end
    
    % Compute the global error
    err(1,index) = norm(xEuler(end)-xExact(end),inf);
    err(2,index) = norm(xRK2(end)-xExact(end),inf);
    
end

figure(1);clf; 
loglog(dt,err,'marker','.','markersize',15,'linestyle','none')
grid on
set(gca,'XDir','reverse')
legend('Euler', 'RK2')
xlabel('dt')
ylabel('Error')