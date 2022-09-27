clc;clear;
close all

%% Initializing the parameters:
tolerance = ... ;
alpha     = ... ;
theta_0   = ... ;

% for plotting purpose:
theta_history = theta_0;
n_plot = 2;

%% Newton iteration:
theta = theta_0;
iteration = true;

while iteration  
    
    % calculate r , dr
    [r,dr] = ...  ;
    
    if abs(r)<tolerance
        break
    else
        % write the direction of dtheta
        dtheta = ... ;
            
        % calculate the new theta     
        theta  = ... ;
            
        % for plotting results
        theta_history(:,n_plot) = theta;
        n_plot = n_plot+1;
    end
end

%% printing the final answer and plotting the function and results:
% Printing final answer in the command window
display(['Final answer is:  '])
fprintf(2,'%f\n',theta)

% plot for 1D:
figure;hold on;
% theta_plot must be changed for each function
theta_plot = -2:0.01:2;
y_plot_main = fungen_r(theta_plot);

plot(theta_plot,y_plot_main);

for i=1:size(theta_history,2)
    y_plot_optim = fungen_r(theta_history(:,i));
    plot(theta_history(:,i),y_plot_optim,'o','color',[1-i/size(theta_history,2) i/size(theta_history,2) 0])
end

