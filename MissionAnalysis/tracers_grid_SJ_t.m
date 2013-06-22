function tr=tracers_grid_SJ_t
n=5;
%% Get the param_struc structure in the figure UserData field.
% fig_file=[folder 'MissionAnalysis/t1.3163_T1.5_points_5.fig'];
fig_file=[folder 'MissionAnalysis/Prove per missione/9luglio/t1.3163_T1.5_5.fig'];
fig=openfig(fig_file);
param_struc=get(fig,'UserData');

% %% get image coordinates
h=get(gca);
% coord_x=get(h.XLabel,'String');
% coord_y=get(h.YLabel,'String');
vx_min=h.XLim(1);
vx_max=h.XLim(2);
vy_min=h.YLim(1);
vy_max=h.YLim(2);

%% Setup the tr structure
tr.n_tracers=n^2;

vx_val=linspace(vx_min,vx_max,n);
vy_val=linspace(vy_min,vy_max,n);
[X,Y]=meshgrid(vx_val,vy_val);
tr.vx=X(:)';
tr.vy=Y(:)';

tr.x=param_struc.x_0*ones(1,n^2);
tr.y=param_struc.y_0*ones(1,n^2);
%tr.e=param_struc.e_0*ones(1,n^2);
%tr.vy=param_struc.vy*ones(1,n^2);
tr.t0=param_struc.t0;
tr.T=param_struc.DT;
tr.mu=param_struc.mu;
tr.ecc=0.04839; % FIXME Attenzione
tr=complete_tracers(tr); % compute nothing