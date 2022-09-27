% ----------------------------------------------------------------------

% In this file you can find the solutions to the Exercises of Matlab
% Tutorial

% ----------------------------------------------------------------------


%% Exercise 2.1:
A = [2,1;-3,1];
disp(det(A))
b = [2;1];
x = A\b;


%% Exercise 2.2:
A = [4.2,0,1.3;1,3.5,pi/2;0,2.2,3];
disp(det(A))
b = [pi;1;2];
x = A\b;


%% Exercise 2.3:
A = [2,0.5,2.5;3,0,3;1,2,3];
disp(det(A)) % so nope


%% Exercise 3.1:
%{ 
% -------------------- copy the following lines (until '%}') to a separate file and save the file:
function a = primeOf(n)
a = [];                % Prepare an empty list.
f = 2;                 % The first possible factor.
while n > 1            % While n still has remaining factors...
    if rem(n,f) == 0   % The remainder of n divided by f might be zero.
        a = [a f];     % If so, it divides n. Add f to the list.
        n = n / f;     % Divide that factor out of n.
    else               % But if f is not a factor of n,
        f = f + 1;     % Add one to f and try again.
    end
end
%Prime factors may be repeated: 12 factors to 2,2,3.
end
%}


%% Exercise 3.2:
fun_gaussian = @(x,mu,sigma)(1/sqrt(2*pi*sigma.^2)*exp(-(x-mu).^2./(2.*sigma.^2)))
x=-6:0.1:6;
mu1=0;           % mean or expectation
sigma1 = sqrt(0.2);      % standard deviation
y1 = fun_gaussian(x,mu1,sigma1);
hold on;         % to allow for multiple plotting
mu2=0;           % mean or expectation
sigma2 = 1;      % standard deviation
y2 = fun_gaussian(x,mu2,sigma2);
mu3=2;           % mean or expectation
sigma3 = sqrt(0.5);      % standard deviation
y3 = fun_gaussian(x,mu3,sigma3);
plot(x,y1)
plot(x,y2)
plot(x,y3)
xlabel('x')
ylabel('p(x)')
legend('\mu=0, \sigma^2=0.2','\mu=0, \sigma^2=1','\mu=2, \sigma^2=0.5')


%% Exercise 4.1:
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


%% Exercise 4.2:
% Define function and successive derivatives
syms x;
f = sin(2*x) + exp(2*x)*cos(x);
Df = jacobian(f,x);
D2f = jacobian(Df,x);
D3f = jacobian(D2f,x);
D4f = jacobian(D3f,x);

% Point x0 around which Taylor approx is done
x0 = 0;
f0 = double(subs(f,x,x0));
Df0 = double(subs(Df,x,x0));
D2f0 = double(subs(D2f,x,x0));
D3f0 = double(subs(D3f,x,x0));
D4f0 = double(subs(D4f,x,x0));

% More terms can be added if needed
f_taylor2 = f0 + Df0*(x-x0) + 1/2*D2f0*(x-x0)^2;
f_taylor3 = f0 + Df0*(x-x0) + 1/2*D2f0*(x-x0)^2 + 1/6*D3f0*(x-x0)^3;
f_taylor4 = f0 + Df0*(x-x0) + 1/2*D2f0*(x-x0)^2 + 1/6*D3f0*(x-x0)^3 + 1/24*D4f0*(x-x0)^4;
fun = matlabFunction(f,'vars',x);
fun_taylor2 = matlabFunction(f_taylor2,'vars',x);
fun_taylor3 = matlabFunction(f_taylor3,'vars',x);
fun_taylor4 = matlabFunction(f_taylor4,'vars',x);

% Plot
clf; hold on;
X = linspace(-2,2,100);
plot(X,fun(X),'k--')
plot(X,fun_taylor2(X))
plot(X,fun_taylor3(X))
plot(X,fun_taylor4(X))

legend('Exact function', 'Taylor order 2', 'Taylor order 3', 'Taylor order 4')

% Notice that the approximation performs better on the left side of the
% plot, whereas the strong nonlinearities of the function makes it
% difficult to approximate on the right side.
