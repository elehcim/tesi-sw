function [integrated_tracers]=integrate_tracers_SOI(tr)
a_jup = 778412027; % km
ecc_jup = 0.04839266;
p=a_jup*(1-ecc_jup^2);
options=odeset('AbsTol',1e-12,'RelTol',1e-12,'Events',@(t,y)(reach_SOI(t,y,tr.mu,tr.ecc,p)));
for i=1:tr.n_tracers
	fprintf('Integrating tracer n. %i\n',i)
	[t,Y]=ode45(...
		@(ni,x)(f_ell(ni,x,tr.mu,tr.ecc)),...
		[tr.t0 tr.t0+tr.T],...
		[tr.x(i); tr.y(i); tr.vx(i); tr.vy(i)],...
		options);
	traj{i,1}=Y;
	traj{i,2}=t;
	% Add energy column
	traj{i,1}(:,5)=0.5*(Y(:,3).^2+Y(:,4).^2)-(Omega(Y(:,1),Y(:,2),tr.mu)./(1+tr.ecc*cos(t)));
end
integrated_tracers=traj;

function [value,isterminal,direction]=reach_SOI(t,y,mu,ecc,p)
%tol = 1e-5;
L=p/(1+ecc*cos(t));
r_SOI=48223000/L; % around 0.0619
value=sqrt((y(1)-(1-mu))^2+y(2)^2)-r_SOI;
isterminal = 1;
direction = 0;