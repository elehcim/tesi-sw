function plot_traj(tr,traj)
for i=1:tr.n_tracers
	hold all
	plot(traj{i,1}(:,1),traj{i,1}(:,2))
	% Initial and final points
	plot(traj{i,1}(1,1),traj{i,1}(1,2),'sg')
	plot(traj{i,1}(end,1),traj{i,1}(end,2),'sr')
	%leg(i)=num2str(i);
end

%% Plot first primary not in scale
plot(-tr.mu,0,'ok','MarkerSize',5)

%% Plot second primary center
plot(1-tr.mu,0,'+k')

% TODO add Legend