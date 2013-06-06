%% Data

jup_grav_par=
a_jup=

t=pi;
e_samples=10;
% energy limits obtained from LCS figures
e_min=
e_max=

e=linspace(e_min,e_max,e_samples);

%% Escape form earth


%% Propagate orbit
t_final=

%% Jupiter injection

dv_jup_inj=deiperbolize(vx,vy,jup_grav_par,t_final,a_jup);

%% plot orbit