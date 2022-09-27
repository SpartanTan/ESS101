clc;clear;

%% Initializing the parameters:
tolerance = ... ;
alpha     = ... ;
theta_0   = ... ;

% for plot
theta_history = theta_0;
n_plot = 2;

%% Newton iteration:
theta = theta_0;
iteration = true;

while iteration         %abs(r) >= tolerance
    
    % calculate r , dr
    [r,dr] = ... ;
    
    if abs(r)<tolerance
        break
    else
        % write the direction of dtheta
        dtheta = ... ;
            
        % calculate the new theta     
        theta = ... ;
            
        % for plot
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
% theta_plot = -12:0.01:5;
theta_plot = (min(min(theta_0),min(theta_history))-1):0.01:(max(max(theta_0),max(theta_history))+1);

for i=1:length(theta_plot)
    y_plot_main = fungen_phi_plot(theta_plot);
end
plot(theta_plot,y_plot_main);

for i=1:size(theta_history,2)
    y_plot_optim = fungen_phi_plot(theta_history(:,i));
    plot(theta_history(:,i),y_plot_optim,'o','color',[1-i/size(theta_history,2) i/size(theta_history,2) 0])
end
