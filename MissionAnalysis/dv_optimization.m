clearvars
addpath('../Tracers/')
%% Data
jup_grav_par = 126711995; % km^3/s^2
sun_grav_par = 132712439935; % km^3/s^2
a_jup = 778412027; % km
ecc_jup = 0.04839266;

%% Select tracers
%tr=select_tracers('../Tracers/Sun_Jupiter_t=220_little.fig');
%load Sun_Jupiter_t=220_little_tracers_20130612-222805.mat
%tr=tracers_grid_SJ_little;

% tr=select_tracers
tr=select_tracers([folder 'MissionAnalysis/Prove per missione/9luglio/t1.3163_T1.5_5.fig'])
% tr=tracers_grid_SJ_t
% load t2.6180_T2_tracers_20130620-131916.mat
%% Escape from earth
%TODO


%% Propagate orbit until Jupiter's SOI
% traj=integrate_tracers_SOI(tr);
load traj_5
tracers_in_SOI=zeros(tr.n_tracers,1);
c=1;
for j=1:tr.n_tracers
	if traj{j,2}(end) ~= tr.t0+tr.T
		tracers_in_SOI(j)=1;
		traj_in_SOI{c,1}=traj{j,1};
		traj_in_SOI{c,2}=traj{j,2};
		c=c+1;
	end
end
n_tracers_in_SOI=sum(tracers_in_SOI);
fprintf('\nNumber of tracers reaching SOI: %i\n', n_tracers_in_SOI);

%%{
%% Plot trajectories
figure
SOI_indexes=find(tracers_in_SOI);
plot_traj(tr.mu,traj_in_SOI,SOI_indexes)
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
nu=zeros(1,tr.n_tracers);
x_syn=zeros(1,tr.n_tracers);
y_syn=zeros(1,tr.n_tracers);
vx_syn=zeros(1,tr.n_tracers);
vy_syn=zeros(1,tr.n_tracers);
dv_jup_inj=nan(1,tr.n_tracers);
fprintf('\n\nJupiter orbit injection delta v:\n')

for j=1:tr.n_tracers
	if tracers_in_SOI(j)
	nu(j)=traj{j,2}(end);
	x_syn(j)=traj{j}(end,1);
	y_syn(j)=traj{j}(end,2);
	vx_syn(j)=traj{j}(end,3);
	vy_syn(j)=traj{j}(end,4);
	[x , y]=r_syn2in(x_syn(j),y_syn(j), nu(j), a_jup, ecc_jup);
	[vx, vy]=v_syn2in(x_syn(j),y_syn(j),vx_syn(j),vy_syn(j),nu(j), ...
		a_jup, ecc_jup, sun_grav_par+jup_grav_par);
	dv_jup_inj(j)=deiperbolize(x,y,vx,vy,jup_grav_par,nu(j),ecc_jup,a_jup);
	fprintf('tracer %02i dv = %.2f km/s\n',j,dv_jup_inj(j))
	end
end

%% choose the min dv and plot that traj
index=find(dv_jup_inj==min(dv_jup_inj));
fprintf('\noptimal tracer n.: %d\ndv = %.2f\n', index,dv_jup_inj(index));
%% Plot trajectories
figure
plot_traj(tr.mu,traj(index),index)
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
%% plot orbit near Jupiter %TODO