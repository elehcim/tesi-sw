function J=cr3bp_jacH(~,y_)
% TODO starting from this code
% mu=0.1;
% x=y_(1);
% y=y_(2);
% Omega_xx= 1-(1-mu)./(((x+mu).^2+y.^2).^(1.5))+(3*(1-mu).*(x+mu).^2)./(((x+mu).^2+y.^2).^(2.5))-mu./(((x-1+mu).^2+y.^2).^(1.5))+(3*mu*(x-1+mu).^2)./(((x-1+mu).^2+y.^2).^(2.5));
% Omega_xy=(3*(1-mu).*(x+mu).*y)./(((x+mu).^2+y.^2).^(2.5))+(3*mu*(x-1+mu).*y)./( ((x-1+mu).^2+y.^2).^(2.5));
% Omega_yy= 1-(1-mu)./(((x+mu).^2+y.^2).^(1.5))+(3*(1-mu).*y.^2)./( ((x+mu).^2+y.^2).^(2.5) )-mu./( ((x-1+mu).^2+y.^2).^(1.5))+(3*mu*y.^2)./( ((x-1+mu).^2+y.^2).^(2.5));
% J = [0, 0, 1, 0;...
% 	0, 0, 0, 1;...
% 	Omega_xx, Omega_xy, 0,2;...
% 	Omega_xy,Omega_yy,-2,0];
end
% Jacobian — Supplying an
% analytical Jacobian often increases the speed and reliability of the
% solution for stiff problems. Set this property to a function FJac,
% where FJac(t,y) computes ∂f/∂y,
% or to the constant value of ∂f/∂y.The Jacobian for the van der Pol Equation (Stiff), described in the MATLAB Mathematics
% documentation, can be coded asfunction J = vdp1000jac(t,y)
% J = [  0                   1
%      (-2000*y(1)*y(2)-1)  (1000*(1-y(1)^2)) ];