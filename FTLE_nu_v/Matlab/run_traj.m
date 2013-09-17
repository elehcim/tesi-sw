%fig_file='../ftle_nuv_mu=0.0010_ecc=0.05_n=200x200_DT=5.000000.fig';
%%
tr=select_tracers;
%%
integrated_tracers=integrate_tracers(tr);
%%
figure
plot_traj(tr.mu,integrated_tracers)