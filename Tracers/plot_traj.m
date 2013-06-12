function plot_traj(tr,traj)
for i=1:tr.n_tracers
	hold all
	plot(traj{i,1}(:,1),traj{i,1}(:,2))
	leg{i}=sprintf('tracer %d',i);
end
legend(leg)

% Initial and final points
for i=1:tr.n_tracers
plot(traj{i,1}(1,1),traj{i,1}(1,2),'sg')
plot(traj{i,1}(end,1),traj{i,1}(end,2),'sr')
end
%% Plot first primary not in scale
plot(-tr.mu,0,'ok','MarkerSize',5)

%% Plot second primary center
plot(1-tr.mu,0,'+k')
axis equal
xlabel('$x$','fontsize',20,'interpreter','latex')
ylabel('$y$','fontsize',20,'interpreter','latex')
set(gca,'FontSize',20)