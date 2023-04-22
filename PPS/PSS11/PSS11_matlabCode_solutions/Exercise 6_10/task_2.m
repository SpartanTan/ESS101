%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%               PSS on Solution to DAE, fully implicit DAE Case
%                   
% Written by Robert Hult, adapted from Sebastien Gros.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear


% WRITE FUNCTIONS ---------------------------------------------------------
subAssignment = 'a';
switch subAssignment
    case 'a'
        z  = sym('z',[1,1],'real');
        x  = sym('x',[1,1],'real');
        dx = sym('dx',[1,1],'real');
        u  = sym('u',[1,1],'real');

        % We know that this DAE has index 1 (very similar to 6.8 in PSS 10)
        F = [dx(1) + tanh(dx(1)) + x(1)*z(1) + u(1);
             tanh(2*u(1)-z(1))-z];
             
        w      = [dx;z]; 
        dF     = jacobian(F,w); 
        x0     = 2;
        tstop  = 10; % termination time of simulation

        un = @(t) 2;
    case 'b' % helicopter dynamics from assignment 1. You can have a look at pages 112-113 in the book to understand the theoretical background.
       
    otherwise
        warning('non-existing subassignment')
        return
end
matlabFunction(F,dF,'File','implicitRHS','vars',{w,x,u});




% SIMULATE ODE ------------------------------------------------------------
tolerance = 1e-6;
printNewtonProgress = true;
z0  = zeros(length(z),1);   % initial guess for newton
dx0 = zeros(length(dx),1); % initial guess for newton
odeWrapper = @(t,x) myImplicitDAEtoODEfunction( t,dx0,x,z0,un(t),tolerance, printNewtonProgress );

tstart  = 0;
odesol  = ode45(odeWrapper,[tstart,tstop],x0);



% PLOT RESULTS ------------------------------------------------------------
% PLOT X
figure(1)
clf
nx= length(x);

if subAssignment == 'a'

    for i = 1:nx
       subplot(nx+1,1,i)
       hold on
       grid on
       plot(odesol.x,odesol.y(i,:))
       xlabel('$t$','interpreter','latex')
       ylabel(sprintf('$x_{%i}$',i),'interpreter','latex')
    end
    % plot |F|
    for i = 1:length(odesol.x)
        t = odesol.x(i);
        x = odesol.y(:,i);
        [~, Ftmp]= myImplicitDAEtoODEfunction( t,dx0,x,z0,un(t),tolerance, false );
        Fnorm(i) = norm(Ftmp);
    end
    subplot(nx+1,1,nx+1)
    semilogy(odesol.x,Fnorm)
    hold on
    grid on

    xlabel('t')
    ylabel('$|F(\dot x,x,z,u)|$','interpreter','latex')

    
elseif subAssignment == 'b'
    
    % plot the positions
    for i = 1:6
        subplot(6,1,i)
        hold on
        grid on
        plot(odesol.x,odesol.y(i,:))
        xlabel('$t$','interpreter','latex')
        ylabel(sprintf('$q_{%i}$',i),'interpreter','latex')
    end
   
    
    figure(2)
    clf;
    
    % plot the speeds
    for i = 1:6
        subplot(6,1,i)
        hold on
        grid on
        plot(odesol.x,odesol.y(6+i,:))
        xlabel('$t$','interpreter','latex')
        ylabel(sprintf('$dq_{%i}$',i),'interpreter','latex')
    end
    
    
    figure(3)
    clf;
    
    % plot |F|
    for i = 1:length(odesol.x)
        t = odesol.x(i);
        x = odesol.y(:,i);
        [~, Ftmp]= myImplicitDAEtoODEfunction( t,dx0,x,z0,un(t),tolerance, false );
        Fnorm(i) = norm(Ftmp);
    end
    semilogy(odesol.x,Fnorm)
    hold on
    grid on

    xlabel('t')
    ylabel('$|F(\dot x,x,z,u)|$','interpreter','latex')

    
else
    
    warning('non-existing subassignment')
    
end
