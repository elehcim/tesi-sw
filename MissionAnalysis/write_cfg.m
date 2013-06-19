clearvars
tempi=[0  0.5236  1.0472  1.5708  2.0944  2.6180  3.1416  3.6652  4.1888  4.7124  5.2360  5.7596];

sun_grav_par = 132712439935; % km^3/s^2
a_jup = 778412027; % km
ecc_jup = 0.04839266;
R = 149600000;	% km
omega= 29.7830/R;	% rad/s

n=300;
T=2;
nvx=n;
nvy=n;
vx_min=[2   2.3  2.3  1.6  0.8  0    0    1.6  1.6  0.8  0   0]
vx_max=[4   4    4    3.6  4    3.3  2.5  4    3.4  4    4   4]


vy_min=[0.5  0    -1.6 -2.6 -3.4 -3.6 -4   -2.5 -2.6 -3.2 -3 -3]
vy_max=[1.4  0.8  -0.4 -1.6 -2.2 -2.4 -2.5 -1.3 -1.5 -2  3 3]


for j=1:length(tempi)
file_name=sprintf('Configuration_%.4f.cfg',tempi(j));

[x(j),y(j)]=earth_synodic(tempi(j), sun_grav_par, a_jup, ecc_jup, omega, R);

data=sprintf(['[parameters]\nmu=9.537e-4	#Sun-Jupiter\necc=0.04839\nDT=%.2f\n',...
	't0=%.4f\n',...
	'n_frames=1\nn_cores=4\nabs_tol=1e-12\nrel_tol=1e-12\n',...
	'file_name=t%.4f\n',...
	'\n[vis.var]\nnvx=%d\n',...
	'vx_min=%.2f\n',...
	'vx_max=%.2f\n',...
	'nvy=%d\n',...
	'vy_min=%.2f\n',...
	'vy_max=%.2f\n',...
	'\n[in.var]\nx=%.4f\ny=%.4f\n',...
	'\n[distance]\nflag=1\nL=778547200\nd1=7e5\nd2=72000\n\n[matlab]\nflag=0'],...
	T,tempi(j),tempi(j),nvx,vx_min(j),vx_max(j),nvy,vy_min(j),vy_max(j),...
	x(j),y(j));

fid=fopen(file_name,'w');
fprintf(fid,data);

end