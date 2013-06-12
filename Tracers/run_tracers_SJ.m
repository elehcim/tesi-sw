clearvars
tr=select_tracers()
% tr=select_tracers('Sun_Jupiter_t=220_little.fig');
% load Sun_Jupiter_t=220_little_tracers_183348  % Tracers con x = x_earth
% load Sun_Jupiter_v_t=220_tracers_171341		% Tracers con impatto
% generate_tracers_SJ
traj=integrate_tracers(tr);
figure
plot_traj(tr,traj)
% size of Jupiter SOF in the adimensionalized system is 0.0619

%% Useful data
% r_sun=6.96e5;		% km
soi_jup=48223000;	% km
L = 778547200;		% km
r_earth_orbit = 149600000;	% km
size_earth_orbit = r_earth_orbit/L;
% size_sun_radius=r_sun/L;
size_jup_soi=soi_jup/L;
% size_jup_sof/size_sun_radius;
%% Plot Jupiter SOI
circle(1-tr.mu,0,size_jup_soi);
%% Plot Earth orbit
circle(-tr.mu,0,size_earth_orbit);