clearvars
% load data_opt_100_grid_dv
% load data_400_grid_dv
load data_2500_grid_dv
%% best 10 dv
% from http://stackoverflow.com/questions/2692482/get-the-indices-of-the-n-largest-elements-in-a-matrix
n=20;
dv_total_copy = dv_total;
for j = 1:n
   [a, Index(j)] = min(dv_total_copy);
   dv_total_copy(Index(j)) = NaN;
end
minimumValues = dv_total(Index)
%%
fi=openfig([folder 'MissionAnalysis/Prove per missione/9luglio/6 to_dv_grid/to_dv_grid_300.fig']);
set(fi,'name','dv scatter')
hold on;
% x_=meshgrid(tr.vx(Index),tr.vy(Index));
% y_=meshgrid(tr.vx(Index),tr.vy(Index));
scatter3(tr.vx(Index),tr.vy(Index),dv_total(Index),80,'r','filled');
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