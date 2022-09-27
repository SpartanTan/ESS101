%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab tips and tricks for ESS101 Modelling and Simulation 2018 
%
% by Robert Hult
%
% The purpose if this script is to demonstrate some of the capabilities of
% Matlab's symbolic toolbox. If you are already familiar with this stuff -
% good for you! :) 
%
% Using these features will greatly simplify your life when dealing with
% the hand-ins of this course.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

% -------------------------------------------------------------------------
%                   DEFINIING SYMBOLIC VARIABLES
% -------------------------------------------------------------------------

% Definition of symbolic vector size [n x 1]
n = 5;
x = sym('x',[n,1]);

% Definition of the same vector with assumption that it is real or positive
x = sym('x',[n,1],'real');
x = sym('x',[n,1],'positive');

% Definition of the same vector with formatted name
x = sym('x_%d',[n,1],'real');

% Definition of symbolic matrix [n x n]
M = sym('m',[n n],'real');

% Definition of symbolic matrix with formatted name
M = sym('m_%d%d',[n n],'real');

% Alternative way of defining symbols.
% Format: symx var1 var2 var3 assumption
syms p g k positive

% -------------------------------------------------------------------------
%                    SYMBOLIC MANIPULATIONS
% -------------------------------------------------------------------------
x = sym('x',[3,1],'real');
syms p k g real
P = [p;k;g];

% You can use matlab functions to build your symbolic expression to the
% desired level of complexity.
f = flipud(x).'*tanh(x);

% you can increase complexity in stages
f = sqrt(f*x(3)-k*x(2)*(cos(x(1))^2+sin(x(1))^2))^k;

% you can simplify expressions
f = simplify(f);

% you can make substitutions 
fsubs = subs(f,k,x(3));         % individual substitution of k with x(3)
fsubs = subs(f,{p,g},{k,x(1)}); % substituting p for k, and g for x(1)
fsubs = subs(f,x,P);            % substitution of vector x with vector P

% you can differentiate w.r.t a specified variable
f_jacobian = jacobian(f,x);
f_gradient = f_jacobian';   % <--- so you don't forget :)

f_hessian  = hessian(f,x);

%% ------------------------------------------------------------------------
%                   FUNCTION GENERATION
% -------------------------------------------------------------------------
% You can generate functions that can be used in matlab from sybolic
% expressions. These can be either anonymous functions or generate matlab
% function files.

% Example:
n = 2;
x = sym('x_%d',[n,1],'real');
x0 = sym('x0_%d',[n,1],'real');
f = (x-x0).^2;

% Anonymous function: gives a handle to f(x,x0).
% NOTE: If yo do not put a vector argument within {}, the function will
% treat the individual vector elements as individual arguments.
f_num = matlabFunction(f,'vars',{x,x0}); 

x0n = rand(2,1);
xn  = rand(2,1);
fn = f_num(xn,x0n);

% You can also write function files (look in the same directory as this
% file for the result)
matlabFunction(f,'file','f_numFile','vars',{x,x0});

% you can call the function like this
fn = f_numFile(xn,x0n);

% You can specify in which folder you want the function to go like this
if not(isdir('functions'))
    mkdir functions % this makes the directory
end
matlabFunction(f,'file','functions/f_numFile','vars',{x,x0});

% you can have multiple outputs to one file. This can be useufull if you,
% for instance, want to get the derivative of a function with the function
% it self
Df  =  jacobian(f,x);
matlabFunction(f,Df,'file','f_numFileWithDerivative','vars',{x,x0});

[fn, Dfn] = f_numFileWithDerivative(xn,x0n);


 
% OBSERVE1: All symbolic variables need to be accounted for when the
% function is generated. That is, the following throws an error that x0_i
% must be included. 
try
f_num = matlabFunction(f,'vars',x); 
catch e
    warning(e.message)
end



% OBSERVE2: In general, the generate functions are much faster than the
% anonymous functions and the way to go!

% OBSERVE3: Using symbolic calculations and code generation removes two of
% the most common implementation errors: Manipulation errors (wrong signs
% etc.) and coding errors (simply writing the wrong thing). 
% Dont be the person that keep on doing these errors after looking through
% this file.


%% ------------------------------------------------------------------------
%                   EXAMPLES
% -------------------------------------------------------------------------
% Some examples on how to use this functionality.

%% Rootfinding with newton.
syms x k real
% we want to find the root of 
f = -atan(x)+atan(2*x) +3*x-pi;
Df = jacobian(f,x);
matlabFunction(f,Df,'file','f_newton','vars',{x});

x = linspace(-5,5);
fx = f_newton(x);

figure(1)
clf
hold on
grid on
plot(x,fx)


xk = 0; % initial guess
it = 0;
fprintf('#it\t|f|\t\t|dx|\n')
while 1
    it = it +1;
    [f,df] = f_newton(xk);
    plot(xk,f,'r*')
    
    xk = xk - f/df;  
    fprintf('%2i\t%1.2e\t%1.2e\n',it,norm(f),norm(f/df))
    if norm(f)<1e-8
        fprintf('Solution is %2.16f\n',xk)
        break
    elseif it > 20
        break
    end
end
 
%% Linearization of a dynamical system xdot = f(x,u)
x = sym('x',[3,1]);
u = sym('u',[2,1]);

x0 = ones(3,1);
u0 = ones(2,1);

f = [-x(1)^2-x(3)*x(2)*tanh(u(1)); -x(2)*x(1)*u(1)*u(2); cos(x(1))*u(1)*sin(u(2))];

DfDx = jacobian(f,x);
DfDu = jacobian(f,u);

% Linearized dynamics
f0 = double(subs(f,[x;u],[x0;u0]));
A  = double(subs(DfDx,[x;u],[x0;u0]));
B  = double(subs(DfDu,[x;u],[x0;u0]));

% NL and Linear systems.
fFun    = matlabFunction(f,'vars',{x,u});
flinFun = @(x,u) f0 + A*(x-x0) + B*(u-u0);


%% ------------------------------------------------------------------------
%                   A PRESENTATION HINT
% -------------------------------------------------------------------------
% When you are writing the report the following function is usefull and
% timesaving

x = sym('x_%d',[2,1],'real');
syms k g real
f = [cos(x(1))*k + sin(x(2))*g; exp(x(1)*x(2))]; 
latexCode = latex(f);

%% ------------------------------------------------------------------------
%                   A WARNING
% -------------------------------------------------------------------------
% While it is very convenient, Matlabs symbolic features can become VERY
% slow for functions of high symbolic complexity. 
return
N  = 3;
dt = 1;

% System
x = sym('x',[3,1]);
f = [10*(x(2)-x(1)); x(1)*(28-x(3))-x(2); x(1)*x(2)-8/3*x(3)];
matlabFunction(f,'file','f_lorentz','vars',{x});

% ERK4 integrator
syms h real

k1 = h*f_lorentz(x);
k2 = h*f_lorentz(x+k1/2);
k3 = h*f_lorentz(x+k2/2);
k4 = h*f_lorentz(x+k3);
xn = x + 1/6*(k1 +2*k2 +2*k3 +k4);


matlabFunction(xn,'file','RK4_lorentz','vars',{x,h});

% N RK4 steps

hnum = dt/N;
xk = x;

% this will take some time to evaluate....
for i = 1:N
    xk = RK4_lorentz(xk,hnum);
end

% try to printout xk(x) if you dare... :)

% Taking the derivative of xk w.r.t x will take some time
Dxk = jacobian(xk,x);

% Writing the functions will take a looong time
matlabFunction(xk,Dxk,'file','N_RK4_lorentz','vars',{x});


% There are other tools that do this in an much more efficient way. 
% A notable mention is "CasADi", which is free of charge and available for
% MATLAB (and several other platforms).


% Good luck coding!

