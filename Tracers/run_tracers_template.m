fig_file='';

%%
tracers=select_tracers(fig_file);
%%
integrated_tracers=integrate_tracers(tracers);
%%
contours=select_contours;
%%
plot_tracers(integrated_tracers,contours)