function [yhat_new,dyhat_new] = yhatFile(in1,u,yhat,in4)
%yhatFile
%    [YHAT_NEW,DYHAT_NEW] = yhatFile(IN1,U,YHAT,IN4)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    28-Sep-2022 14:56:30

a = in1(1,:);
b = in1(2,:);
dyhat1 = in4(:,1);
dyhat2 = in4(:,2);
t2 = -yhat;
yhat_new = a.*t2+b.*u;
if nargout > 1
    dyhat_new = [t2-a.*dyhat1,u-a.*dyhat2];
end
