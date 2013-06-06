function dy = dg( t,y )
%DG Double Gyre dynamical system
%  x_dot = dg( t,y )

epsilon=0.25;  %% amount of sideways perturbation of gyre %%
omega=2*pi/10;
a=epsilon*sin(omega*t);
b=1-2*epsilon*sin(omega*t);

A=0.1;   %% Amplitude of streamfunction %%

f=a*y(1)^2+b*y(1);
u=-pi*A*(sin(pi*f))*cos(pi*y(2));
v=pi*A*cos(pi*f)*sin(pi*y(2))*(2*a*y(1)+b);

dy = [u;v];
end