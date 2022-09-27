clear all
clc

%% Simulator ARX as model predictor

syms a b yhat dyhat u real

% two parameters
theta = [a;b];
% we also define symbolically dyhat_dtheta
dyhat = sym('dyhat',[1,2]);

% We define the ARX simulator at a generic time (new), function of the
% simulated output at the previous time (yhat).
yhat_new  = -a*yhat + b*u;

% Now we need to define the derivative of the predictor yhat_new wrt the
% parameters (note that yhat_new contains yhat, which is also function of the parameters,
% but we defined symbolically this derivative, that we can use here).
dyhat_new = ...  ;

% now we export the function providing the simulator output and its
% sensitivity wrt the parameters.
matlabFunction(yhat_new,dyhat_new, 'file','yhatFile','vars',{theta,u,yhat,dyhat});