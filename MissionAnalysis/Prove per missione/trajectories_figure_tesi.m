addpath(['../../Tracers/'])
addpath(['../'])
tr=select_tracers('9luglio/')
traj=integrate_tracers(tr);
figure
plot_traj(tr.mu,traj)