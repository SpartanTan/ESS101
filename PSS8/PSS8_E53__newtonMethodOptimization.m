clear all
close all
clc

Order_u = 1; % number of b parameters
Order_y = 1; % number of a parameters

% define number of data, Newton tolerance and alpha
Ndata = 1e3;
tol   = 1e-6;
alpha = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% DATA GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% True value of the parameters
theta.true = [0.95;1];

% Generate data
y.data = zeros(max(Order_u,Order_y)+Ndata,1);
u      = random('norm',0,1,  [max(Order_u,Order_y)+Ndata,1]); % random noise with std=1
e      = random('norm',0,0.1,[max(Order_u,Order_y)+Ndata,1]); % random noise with std=0.1

% Simulate system to obtain the DATA
for k = 2:Ndata+1
    y.data(k) = yhatFile(theta.true,u(k),y.data(k-1),zeros(1,2)) + e(k);
end

figure(1);clf;hold on
plot(y.data(Order_y+1:end),'marker','.','linestyle','none','color','b','markersize',20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PARAMETER  ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the fitting function
theta.est = [0.5;1]; % initial condition for Newton's iterations
iter      = true;
y.hat     = 0;
niter = 1;
while iter
    % initialize variables for Jacobians and Hessian
    dJ = zeros(1,2);
    d2J = zeros(2,2);
    dyhat = zeros(1,2);
    
    % Simulate model for current theta.est
    for k = 2:Ndata+1
        [y.hat(k),dyhat] = yhatFile(theta.est,u(k),y.data(k-1),dyhat);
        % the above function provides the simulated output y.hat and the
        % sensitivity dyhat for 1 time instant. Hence we can build jacobian
        % and Hessian of the cost function for all time instants (summing all together)
        % Jacobian :
        dJ  = dJ + ...  ;
        % Hessian (here we are going to use an approximation):
        d2J = d2J + ...  ;
    end
    
    % we print the simulator for the current value of the parameters,
    % comparing it with the data.
    figure(1);clf;hold on
    plot(y.data(Order_y+1:end),'marker','.','linestyle','none','color','b','markersize',20)
    plot( y.hat(Order_y+1:end),'marker','.','linestyle','none','color','r','markersize',20)
    drawnow
    
    % Here we have to take our Newton step
    dtheta        = ...  ;
    theta.est     = theta.est + alpha*dtheta; 
    
    % Can we stop iterating?
    ExitCriterion(niter) = norm(dJ,inf);
    if ExitCriterion(end) < tol 
        iter = false;
    else
        niter = niter + 1;
    end
   
end

% Plot convergence of the Newton iteration
figure(2)
semilogy(ExitCriterion,'marker','.','linestyle','none','color','k','markersize',20)
grid on

fprintf('true theta = \n        %.4f , %.4f \n',theta.true(1),theta.true(2));
fprintf(2,'estimated theta = \n        %.4f , %.4f \n',theta.est(1),theta.est(2));