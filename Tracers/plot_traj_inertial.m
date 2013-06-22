function plot_traj_inertial( traj, index)
a = 778412027; % km
r_earth_orbit = 149600000;	% km
ecc = 0.04839266;

hold all
ellipse(a,a*sqrt(1-ecc^2),0,0,0)
circle(0,0,r_earth_orbit)
plot(0,0,'+')

n=size(traj,1);
for i=1:n
	x_syn=traj{i}(:,1);
	y_syn=traj{i}(:,2);
	nu=traj{i,2};
	[x_in,y_in]=r_syn2in(x_syn,y_syn,nu, a,ecc);
	plot(x_in,y_in)
	leg{i}=sprintf(' %d',i);
end
if nargin == 3
	for j=1:length(index) % Useful to plot the chosen value when optimizing
	leg{j}=sprintf(' %d',index(j));
	legend(leg,'fontsize',10)
	end
else
	legend(leg,'fontsize',10)
end
xlabel('$x$','fontsize',20,'interpreter','latex')
ylabel('$y$','fontsize',20,'interpreter','latex')
set(gca,'FontSize',20)
axis equal
end

