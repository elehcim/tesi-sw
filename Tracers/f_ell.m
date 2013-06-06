function x_out=f_ell(ni,x,mu,e)
K=1+e*cos(ni);
x_prime=x(3);
y_prime=x(4);
wx_prime=2*x(4)+(x(1)-((1-mu).*(x(1)+mu))./(((x(1)+mu).^2+x(2).^2).^(1.5))-(mu.*(x(1)-1+mu))./(((x(1)-1+mu).^2+x(2).^2).^(1.5)))/K;
wy_prime=-2*x(3)+(x(2)-(1-mu).*x(2)./(((x(1)+mu).^2+x(2).^2).^(1.5))-(mu*x(2))./(((x(1)-1+mu).^2+x(2).^2).^(1.5)))/K;
x_out=[x_prime, y_prime, wx_prime, wy_prime]';