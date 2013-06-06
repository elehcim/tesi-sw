function om = Omega_mas( x,y,mu )
% Omega_mas( x , y , mu )
% Omega_mas Compute the negative effective potential following Masdemont
% article
r1=sqrt((x-mu).^2+y.^2);
r2=sqrt((x+1-mu).^2+y.^2);
om = 0.5.*(x.^2+y.^2) + (1-mu)./r1 + mu./r2 + 0.5.*(mu.*(1-mu));
end