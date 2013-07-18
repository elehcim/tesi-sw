tr=select_tracers([folder '../../../../Desktop/Figure_tesi/'])
traj=integrate_tracers(tr);
figure
plot_traj(tr.mu,traj)