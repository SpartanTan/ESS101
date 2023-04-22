clear all
clc

%% OE as model predictor

syms a b yhat dyhat u real

% two parameters
theta = [a;b];
% we also define symbolically dyhat_dtheta
dyhat = sym('dyhat',[1,2]);

% We define the OE at a generic time (new), function of the
% simulated output at the previous time (yhat).
yhat_new  = -a*yhat + b*u;

% Now we need to define the derivative of the predictor yhat_new wrt the
% parameters (note that yhat_new contains yhat, which is also function of the parameters,
% but we defined symbolically this derivative, that we can use here).
dyhat_new = jacobian(yhat_new, theta) +  jacobian(yhat_new, yhat)*dyhat;  %(5.9)

% now we export the function providing the simulator output and its
% sensitivity wrt the parameters.
matlabFunction(yhat_new,dyhat_new, 'file','yhatFile','vars',{theta,u,yhat,dyhat});