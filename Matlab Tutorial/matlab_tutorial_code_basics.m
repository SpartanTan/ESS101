%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab tips and tricks for ESS101 Modelling and Simulation 2019 
%
% by Ramin and Rémi
%
% > The purpose of this script is to introduce some of the most basic and 
% useful capabilities of Matlab. If you are already familiar with this 
% stuff - good for you! :) 
%
% > Please have a look at the pdf file that comes with this script for
% further info. 
% 
% > And please note that this script is not meant to be executed as a 
% whole, but is rather a convenient way to gather useful functions that
% you may want to use during the course.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear


% -------------------------------------------------------------------------
%                   VECTORS AND MATRICES
% -------------------------------------------------------------------------

% Use the help function whenever you want to have more info on what
% something does
help arctan;

% You may alternatively use the doc function
doc arctan;

% It is easy to print anything with disp. It works with most variable types
disp('hello, world')

% Define vectors and matrices by hand
a = [42, 10, 7.5, 88, 3.14];                        % line vector
b = [42 10 7.5 88 3.14];                            % equivalent to a
c = [42; 10; 7.5; 88; 3.14];                        % column vector
A = [1, 2, 3; 4, 5, 6];                             % 2*3 matrix
B = [0.42, 0.42; 0.28, 0.28; 0.1, 0.1; 0.31, 0.29]; % 4*2 matrix

% Check the dimensions of a vector/matrix
size(a)
length(a) % only use length for vectors !
size(B)
length(B) % see?

% Define a sequences of numbers as arrays (i.e. vectors)
a = 0:2:20;          % all even integers between 0 and 20
b = 0:0.1:1; 
c = linspace(0,1,11) % has 11 values. often used for plotting figures

% Generate vectors/matrices containing only zeros or ones (or any scalar)
a = zeros(10,1);
b = ones(10,1);
A = zeros(4,2);
B = 41*ones(3,3);
I = eye(10);      % generate identity matrix
 
% + and - operations between vectors and matrices are written in the same 
% way as with scalars 
b = b - 2;
c = a - b + 2;
B = B + 1;
C = 2*B - (B + B);
C = A + B;         % matrices being added should have the same dimension !

% You can access single matrix elements A(i,j) and vector elements a(i) 
% just as you would imagine
A = [1, 2, 3; 4, 5, 6];
a = 1:10;
A(1,3)
a(5)

% You can even select only some parts of a matrix/vector
B = A(:,1);   % select first column of A
B = A(2,:);   % select second row of A
B = A(1,2:3); % elements on first row of A and belonging to columns 2 or 3
b = a(2:5);
b = a(4:end); % all elements from element 4 to the end of a

% Transposition works like this
trans_A = A';
trans_a = a';


% -------------------------------------------------------------------------
%                   LINEAR ALGEBRA
% -------------------------------------------------------------------------

% The multiplication of matrices is done with the operator '*'
A = ones(2,3);
B = rand(3,3);    % generate a matrix with random elements between 0 and 1
C = A*B;
D = B*A;          % mind the dimensions when using matrix multiplication !
E = pi*B;         % but no problem for multiplying with a scalar

% Please note that the previous is different from term-by-term
% multiplication, which can be done with operator '.*'
I = eye(3);
C = B*I;
D = B.*I; % notice the difference with C
E = I./B; % also works for term-by-term division

% Use 'det' to get the determinant of any (square) matrix
A = [1, 3; 4, 2];
d = det(A);

% If (AND ONLY IF) the determinant is not 0 (or very close to 0) may you
% use 'inv' to get the inverse of the matrix
inv_A = inv(A);

% You can even play with the eigenvalues of the matrix if you want
lambda = eig(A);

% Use operator '\' to solve linear systems of the form Ax = b. But first
% make sure that A is invertible. Don't be that student that tries to
% invert non-invertible matrices ;)
A = [1 4 5; 2 3 5; 5 9 8];
b = [1; 1; 1];
x = A\b;


% -------------------------------------------------------------------------
%                   BUILT-IN FUNCTIONS
% -------------------------------------------------------------------------

help elfun   % get a list of elementary functions
help specfun % get a list of special functions
 
% Trigonometric and exponential functions are built-in just as you would 
% expect, no surprises here. Same for pi and i (but not e)
cos(pi);
sin(pi/2);
tan(pi/4);
exp(1);
log(exp(1));
sqrt(2);
2^8;         % x to the power n is written x^n 

% Some other useful functions whose names give them away
a = [16; -23; 10; 42; -37; 22; -30];
sum(a);
[x,i] = min(a); % x is the minimum value and i is its index in a
[x,i] = max(a); 
a = abs(a);     % works with matrices, vectors and scalars
round(3.66);    % round to nearest integer
sign(-2);

% -------------------------------------------------------------------------
%                   CONTROL FLOW
% -------------------------------------------------------------------------
% if:
a = rand(1,1);
if a>0.5
    a=1;
else
    a=0;
end

% ---------------------------------------
% switch, example 1
a = 3;
switch a
    case 1
        b = 1;
    case 2
        b = 10;
    case 3
        b = 100;
    otherwise
        b = 0;
end
% switch, example 2
a = 'simulate';
switch a
    case 'simulate'
        disp('simulation is running ...');
    case 'predict'
        disp('prediction is running ...');
    otherwise
        disp('error!');
end

% ---------------------------------------
% for  ->   to create loops
x = rand(1,10);     % create a vector containing rand numbers
y = zeros(1,10);
for i = 1:length(x)
    if x(i)>0.5
        y(i) = 1;
    else
        y(i) = -1;
    end
end

% ---------------------------------------
% while -> to create loops
time = 0;
simulationRun = true;
while simulationRun
    time = time + 0.1;
    if time >= 10
        break
    end
end

% -------------------------------------------------------------------------
%                   DATA IMPORT AND EXPORT
% -------------------------------------------------------------------------
x=1;y=2;
save filename           % Saves all workspace variables to 'filename.mat'
save filename x,y       % Saves variables 'x' and 'y' in 'filename.mat'
% ---------------------------------------
load filename           % Loads all variables from the file 'filename'
load x                  % Loads only the variable 'x' from the file

% -------------------------------------------------------------------------
%                   2D PLOTTING
% -------------------------------------------------------------------------
% plot sin(x) funtion for x between -2pi and pi:
x=-pi:0.1:pi;
y1=sin(x);
plot(x,y1)
title(['Example of plot function'])
xlabel('xlabel text')
ylabel('ylabel text')
hold on
grid
y2 = cos(x);
plot(x,y2,'r--','LineWidth',2)
legend('sin(x)','cos(x)')

% -------------------------------------------------------------------------
%                   FUNCTION DEFINITION
% ------------------------------------------------------------------------- 
% 1) inline function definition:
functionHandle = @(input)(input*2+sin(input));
functionHandle(3);

% 2) defining function in a separate mFile:
% please copy the following code in a new mFile and save it with the
% function Name:
%{
function [c1,c2] = functionHandle(a1,a2)
c1 = a1+a2;
c2 = a1-a2;
end
%}

% 3) Using the symbolic toolbox funtion "matlabFunction"


