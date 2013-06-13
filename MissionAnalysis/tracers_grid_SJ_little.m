function tr=tracers_grid_SJ_little
%% Get the param_struc structure in the figure UserData field.
fig_file='../Tracers/Sun_Jupiter_t=220_little.fig';
fig=openfig(fig_file);
param_struc=get(fig,'UserData');

% %% get image coordinates
h=get(gca);
% coord_x=get(h.XLabel,'String');
% coord_y=get(h.YLabel,'String');
x_min=h.XLim(1);
x_max=h.XLim(2);
vx_min=h.YLim(1);
vx_max=h.YLim(2);

%% Setup the tr structure
n=5;
tr.n_tracers=n;
tr.x=linspace(x_min,x_max,tr.n_tracers);
%tr.x=-0.1922*ones(1,n);
tr.y=param_struc.y_0*ones(1,n);
tr.vx=linspace(vx_min,vx_max,tr.n_tracers);
%tr.vx=linspace(1.534,1.789,n); % these value are proper of the LCS
tr.e=param_struc.e_0*ones(1,n);
tr.t0=param_struc.t0;
tr.T=param_struc.DT;
tr.mu=param_struc.mu;
tr.ecc=0.04839; % FIXME Attenzione
tr=complete_tracers(tr); % compute vy
disp(tr)