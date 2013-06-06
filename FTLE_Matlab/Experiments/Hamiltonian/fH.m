function x_dot = fH(~,x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mu=0.1;
x_dot=zeros(4,1);
x_dot(1,1)=x(3)+x(2);
x_dot(2,1)=x(4)-x(1);
x_dot(3,1)=x(4)-x(1)+x(1)-((1-mu).*(x(1)+mu))./(((x(1)+mu).^2+x(2).^2).^(1.5))-(mu.*(x(1)-1+mu))./(((x(1)-1+mu).^2+x(2).^2).^(1.5));
x_dot(4,1)=-x(3)-x(2)+x(2)-(1-mu).*x(2)./(((x(1)+mu).^2+x(2).^2).^(1.5))-(mu*x(2))./(((x(1)-1+mu).^2+x(2).^2).^(1.5));
end

