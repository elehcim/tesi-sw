clearvars
addpath('../Plot_ftle/')
% load data_100_grid_dv
% load data_400_grid_dv
load data_2500_grid_dv
%% best n dv
% from http://stackoverflow.com/questions/2692482/get-the-indices-of-the-n-largest-elements-in-a-matrix
n=50;
dv_total_copy = dv_total;
for j = 1:n
   [a, Index(j)] = min(dv_total_copy);
   dv_total_copy(Index(j)) = NaN;
end
minimumValues = dv_total(Index);

%% Isolate times
for j=1:tr.n_tracers
transfer_time(j)=traj{j,2}(end);
end
for j=1:n_tracers_in_SOI
transfer_time_in_SOI(j)=traj_in_SOI{j,2}(end);
end

%% correlation dv vs. transfer time
DV=dv_total(logical(tracers_in_SOI));
T=transfer_time_in_SOI;
scatter(T,DV)
xlabel('Transfer time (months)','interpreter','latex','fontsize',35)
ylabel('$\Delta v$ (km/s)','interpreter','latex','fontsize',35)
%% correlation between velocity vector and dv
for j=1:n_tracers_in_SOI
Velocity(j)=sqrt(traj{j,1}(1,3)^2+(traj{j,1}(1,4))^2);
end
scatter(Velocity,DV)
xlabel('$|V|$ (km/s)','interpreter','latex','fontsize',35)
ylabel('$\Delta v$ (km/s)','interpreter','latex','fontsize',35)
%% All reaching points
fi=openfig([folder 'MissionAnalysis/Prove per missione/9luglio/6 to_dv_grid/to_dv_grid_300.fig']);
set(fi,'name','dv span')
hold on;
% x_=meshgrid(tr.vx(Index),tr.vy(Index));
% y_=meshgrid(tr.vx(Index),tr.vy(Index));
color_span=(dv_total-min(dv_total))/(max(dv_total)-min(dv_total))*(2.2-1.1)+1.1;
scatter(tr.vx,tr.vy,25,...
	color_span,'filled');
grid on
hold off
%% only best n with lines
fi=openfig([folder 'MissionAnalysis/Prove per missione/9luglio/6 to_dv_grid/to_dv_grid_300.fig']);
set(fi,'name','dv scatter')
hold on;
% x_=meshgrid(tr.vx(Index),tr.vy(Index));
% y_=meshgrid(tr.vx(Index),tr.vy(Index));
scatter3(tr.vx(Index),tr.vy(Index),dv_total(Index),80,...
	dv_total(Index),'filled');
zlabel('$\Delta v$','interpreter','latex','fontsize',35)
for j=1:n
plot3([tr.vx(Index(j)),tr.vx(Index(j))],...
	[tr.vy(Index(j)),tr.vy(Index(j))],...
	[0,dv_total(Index(j))],'k', 'linewidth',3)
end
box off
grid on
pbaspect([1,1,.2])
hold off

%% all with lines
fi=openfig([folder 'MissionAnalysis/Prove per missione/9luglio/6 to_dv_grid/to_dv_grid_300.fig']);
set(fi,'name','dv scatter 3d')
hold on;
% x_=meshgrid(tr.vx(Index),tr.vy(Index));
% y_=meshgrid(tr.vx(Index),tr.vy(Index));
color_span=(dv_total-min(dv_total))/(max(dv_total)-min(dv_total))*(2.2-1.1)+1.1;
scatter3(tr.vx,tr.vy,color_span, 80,'r','filled');
zlabel('$\Delta v$','interpreter','latex','fontsize',35)
for j=1:tr.n_tracers
plot3([tr.vx(j),tr.vx(j)],...
	[tr.vy(j),tr.vy(j)],...
	[0,color_span],'k', 'linewidth',3)
end
% FIXME adjust z axis
box off
grid on
pbaspect([1,1,.2])
hold off
%% Histogram
figure('name','dv total')
bar(dv_total(Index))
set(gca,'XTickLabel',Index)

%%
%%
% figure('name','tracers'' fingerprint')
% fi1=openfig([folder 'MissionAnalysis/Prove per missione/9luglio/6 to_dv_grid/to_dv_grid_300.fig']);
% plot(gca,tr.vx,tr.vy,'+k','MarkerSize',2);
% % for j=1:tr.n_tracers
% % text(tr.vx,tr.vy,sprintf(' %d',j),'fontsize',3,'VerticalAlignment','top')
% % end