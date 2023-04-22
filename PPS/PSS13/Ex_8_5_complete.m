%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE BUTCHER TABLEU FOR GAUSS-LEGENDRE METHOD (EXERCISE 8.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
s = 4; % number of stages
fprintf('1) stages are selected to s = %i\n\n',s)


% 2) ----------------------------------------------------------------------
fprintf('2) Find the grid points tau as the roots of the Legendre polynomial:\n\nPs = \n')
syms tau real
P  = (tau^2-tau)^s;
Ps = (1/factorial(s))*diff(P,s);
pretty(Ps)

rootsOfPs = sort(solve(Ps));
fprintf('These are:\n\n')
for i =1:length(rootsOfPs)
    fprintf('r%i \t = ',i) 
    disp(rootsOfPs(i))
end
collocationTimes = rootsOfPs;


% 3) ----------------------------------------------------------------------
fprintf('3) Build Lagrange Polynomials:\n \n')

for i = 1:s
    l(i,1) = sym(1);
    for j = 1:s
        if not(j==i)
            l(i,1) = l(i,1)*(tau-collocationTimes(j))/(collocationTimes(i)-collocationTimes(j));
        end
    end
end


% 4) ----------------------------------------------------------------------
fprintf('4) Build Integrals of Lagrange Polynomial :\n \n')
for i = 1:length(l)
    L(i) = int(l(i),[0,tau]);
end


% 5) ----------------------------------------------------------------------
fprintf('5) Evaluate integrals to build A,b,c :\n \n')

% FOR A
A = sym(zeros(s,s));
for i = 1:s
    for j=1:s
        A(i,j) = simplify(subs(L(j),tau,collocationTimes(i)));
    end
end

% FOR b
b = sym(zeros(1,s));
for i = 1:s
    b(i) = simplify(subs(L(i),tau,1));
end

% FOR c
c = sym(zeros(s,1));
for i = 1:s
   c(i) = simplify(collocationTimes(i));
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
