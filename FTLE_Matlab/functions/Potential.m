function Pot = Potential(x,y,mu)
% Computation of effective potential (Short pag. 15)
r1=sqrt((x+mu).^2+y.^2);
r2=sqrt((x-1+mu).^2+y.^2);
Pot=(1-mu)./r1+mu./r2+0.5*(x.^2+y.^2);
end

