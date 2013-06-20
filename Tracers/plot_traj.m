function plot_traj(mu,traj,index)
if nargin < 3
	index = 0;
end
a=size(traj);
n=a(1);
for i=1:n
	hold all
	plot(traj{i,1}(:,1),traj{i,1}(:,2))
	leg{i}=sprintf(' %d',i);
end
if nargin == 3 % Useful to plot the chose value when optimizing
	single_leg=sprintf(' %d',index);
	legend({single_leg},'fontsize',10)
else
	legend(leg,'fontsize',10)
end

% Initial and final points
for i=1:n
plot(traj{i,1}(1,1),traj{i,1}(1,2),'sg')
plot(traj{i,1}(end,1),traj{i,1}(end,2),'sr')
end
%% Plot first primary not in scale
plot(-mu,0,'ok','MarkerSize',5)

%% Plot second primary center
plot(1-mu,0,'+k')
axis equal
xlabel('$x$','fontsize',20,'interpreter','latex')
ylabel('$y$','fontsize',20,'interpreter','latex')
set(gca,'FontSize',20)