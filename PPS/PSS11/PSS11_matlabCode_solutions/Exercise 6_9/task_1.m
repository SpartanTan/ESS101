%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%               PSS on Solution to DAE, semi-explict DAE Case
%                   
% Written by Robert Hult, adapted from Sebastien Gros.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

% WRITE FUNCTIONS ---------------------------------------------------------
subAssignment = 'd';
switch subAssignment
    case 'b'
        z = sym('z',[2,1],'real');
        x = sym('x',[2,1],'real');
        u = sym('u',[1,1],'real');

        f = [z(1)^2+x(2); 
             z(2)-x(1)];
        g = [z(1)^3+z(1)+x(2)+1;
             z(2)^3+z(2)-x(1)-tanh(z(1))];
        
        Dg     = jacobian(g,z); 
        x0     = [4; 2];
        tstop  =  10; % termination time of simulation

        un = @(t) 0;    % this only gives a function of t which is zero
    case 'c'
        z = sym('z',[1,1],'real');
        x = sym('x',[2,1],'real');
        u = sym('u',[1,1],'real');
        
        f = [z*x(1);
             z*x(2)-1];
        g      = x(1)^2+x(2)^3-1;
        Dg     = jacobian(g,z);
        x0     = [1; 0];
        
        useReducedForm = true; % false to get the first problem
        if useReducedForm 
            useAlternateInitialcondition = true;
            if useAlternateInitialcondition % choosing another IC allows the solution to hit t= 10 before det(jacobian(g,z)) = 0 (this doesn't remove the problem, it only moves it to later)
                x0    = [0; 0.999999999];
                x0(1) = min(double(solve(subs(g,x(2),x0(2)),x(1)))); % to get a consistent IC (consistency condition)
            end
            g  = 2*x(1)*f(1)+3*x(2)^2*f(2); % reducing the index of the previous DAE system by 1. but we still have not reached index 1.     
            Dg = jacobian(g,z); 
        end
        tstop  =  30; % termination time of simulation
        
        un = @(t) 0;
        

    case 'd' 
        % The issue here is the same as the previous assignment (look at the determinant of dg in the printout)
        
        z = sym('z',[2,1],'real');
        x = sym('x',[3,1],'real');
        u = sym('u',[1,1],'real');

        f = [  10*x(1) * z(2);
            -30*z(2)*x(2)^2 - z(1)^2 - 10;
            z(1) - 5*z(2)];
        
        g = [x(3) + z(1) - 5*z(2);
            5*z(2) - z(1) - x(3) + 6*x(2)^2 *z(1)^2 + 20*x(1)^2*z(2) + 180*x(2)^4*z(2) + x(1)^2 + 60*x(2)^2-2*x(2)^3 ];
        
        Dg     = jacobian(g,z);
        x0     = [4; 2; 0];
        tstop  =  0.6; % termination time of simulation
        

        un = @(t) 0;
    otherwise
        warning('non-existing subassignment')
end
matlabFunction(g,Dg,'File','algebraicRhs','vars',{x,z,u});
matlabFunction(f,'File','differentialRhs','vars',{x,z,u});


% SIMULATE ODE ------------------------------------------------------------
tolerance = 1e-8;
printNewtonProgress = true;
z0 = zeros(length(z),1); % initial guess for newton
odeWrapper = @(t,x) myDAEtoODEfunction( t,x,z0,un(t),tolerance, printNewtonProgress );

tstart = 0;
odesol  = ode45(odeWrapper,[tstart,tstop],x0,odeset('RelTol',1e-8,'AbsTol',1e-8));



% PLOT RESULTS ------------------------------------------------------------
% PLOT X
figure(1)
clf
nx= length(x);

for i = 1:nx
   subplot(nx+2,1,i)
   hold on
   grid on
   plot(odesol.x,odesol.y(i,:))
   xlabel('t')
   ylabel(sprintf('x_%i',i))
end
% plot |g|
for i = 1:length(odesol.x)
    t = odesol.x(i);
    x = odesol.y(:,i);
    [~, gtmp,z,detdg(i)]= myDAEtoODEfunction( t, x, z0,un(t),tolerance, false );
    
    gnorm(i) = norm(gtmp); 
end
subplot(nx+2,1,nx+1)
semilogy(odesol.x,gnorm)
hold on
grid on
xlabel('t')
ylabel('$|g(x,z)|$','interpreter','latex')

subplot(nx+2,1,nx+2)
semilogy(odesol.x,abs(detdg))
hold on
grid on
xlabel('t')
ylabel('$\mathrm{det}\left(\frac{\partial g }{\partial z}\right)$','interpreter','latex')




