clearvars
addpath('../Tracers/')
%% Data
jup_grav_par = 126711995; % km^3/s^2
a_jup = 778412027; % km
ecc_jup = 0.04839266;

%% Select tracers
tr=select_tracers('../Tracers/Sun_Jupiter_t=220_little.fig');
%% Escape form earth
%TODO
%load Sun_Jupiter_t=220_little_tracers_20130612-222805.mat

%% Propagate orbit until Jupiter's SOI
traj=integrate_tracers_SOI(tr);
%{
%% Plot trajectories
figure
plot_traj(tr,traj)
% size of Jupiter SOF in the adimensionalized system is 0.0619
% Useful data
% r_sun=6.96e5;		% km
soi_jup=48223000;	% km
L = 778547200;		% km
r_earth_orbit = 149600000;	% km
size_earth_orbit = r_earth_orbit/L;
% size_sun_radius=r_sun/L;
size_jup_soi=soi_jup/L;
% size_jup_sof/size_sun_radius;
% Plot Jupiter SOI
circle(1-tr.mu,0,size_jup_soi);
% Plot Earth orbit
circle(-tr.mu,0,size_earth_orbit);
%}
%% Jupiter injection
for j=1:tr.n_tracers
	nu(j)=traj{j,2}(end);
	[x , y]=r_syn2in(traj{j}(end,1),traj{j}(end,2), nu(j), a_jup, ecc_jup);
	[vx, vy]=v_syn2in(traj{j}(end,1),traj{j}(end,2),...
		traj{j}(end,3),traj{j}(end,4),nu(j), a_jup, ecc_jup, jup_grav_par);
	dv_jup_inj(j)=deiperbolize(x,y,vx,vy,jup_grav_par,nu(j),ecc_jup,a_jup);
end
%% plot orbit