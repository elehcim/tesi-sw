function om = Omega( x,y,mu )
% Omega( x , y , mu )
% Omega Compute the negative effective potential
r1=sqrt((x+mu).^2+y.^2);
r2=sqrt((x-1+mu).^2+y.^2);
om = 0.5.*(x.^2+y.^2)+(1-mu)./r1+mu./r2+0.5.*(mu.*(1-mu));
end