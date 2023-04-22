%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE BUTCHER TABLEU FOR GAUSS-LEGENDRE METHOD (EXERCISE 8.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
s = 2; % number of stages
fprintf('1) stages are selected to s = %i\n\n',s)


% 2) ----------------------------------------------------------------------
fprintf('2) Find the grid points tau as the roots of the Legendre polynomial:\n\nPs = \n')
syms tau real
% Ps = .... < ---------------- YOUR STUFF HERE
% pretty(Ps)

% rootsOfPs =  ...
fprintf('These are:\n\n')
for i =1:length(rootsOfPs)
    fprintf('r%i \t = ',i) 
    disp(rootsOfPs(i))
end
collocationTimes = rootsOfPs;


% 3) ----------------------------------------------------------------------
fprintf('3) Build Lagrange Polynomials:\n \n')

for i = 1:s
    % YOUR STUFF HERE 
end


% 4) ----------------------------------------------------------------------
fprintf('4) Build Integrals of Lagrange Polynomial :\n \n')
for i = 1:s
%     L(i) = .......
end


% 5) ----------------------------------------------------------------------
fprintf('5) Evaluate integrals to build A,b,c :\n \n')

% FOR A
A = sym(zeros(s,s));
for i = 1:s
    for j=1:s
%         A(i,j) = ... <-------------- YOUR STUFF HERE
    end
end

% FOR b
b = sym(zeros(1,s));
for i = 1:s
%     b(i) = ... <-------------- YOUR STUFF HERE
end

% FOR c
c = sym(zeros(s,1));
for i = 1:s
%    c(i) = ... <-------------- YOUR STUFF HERE
end

A
b
c


% SHOULD BE FOR s = 2
% A =
%  
% [             1/4, 1/4 - 3^(1/2)/6]
% [ 3^(1/2)/6 + 1/4,             1/4]
%  
%  
% b =
%  
% [ 1/2, 1/2]
%  
%  
% c =
%  
%  1/2 - 3^(1/2)/6
%  3^(1/2)/6 + 1/2

butcher.A = eval(A);
butcher.b = eval(b);
butcher.c = eval(c);
save('myButcherTableu.mat','butcher')
