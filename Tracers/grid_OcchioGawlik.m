% generate grid of points for Gawlik
clearvars
n=10;
tr.n_tracers=n^2;

fig_file='OcchioGawlik.fig';
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

x_val=linspace(x_min,x_max,n);
vx_val=linspace(vx_min,vx_max,n);
[X,Y]=meshgrid(x_val,vx_val);
tr.x=X(:)';
tr.vx=Y(:)';

% tr.x=linspace(-0.44,-0.15,n);
tr.y=zeros(1,n^2);
% tr.vx=linspace(0.44,1.75,n);
tr.e=-1.74*ones(1,n^2);
tr.t0=0;
tr.T=2;
tr.mu=0.1;
tr.ecc=0.04;
tr=complete_tracers(tr); % compute vy
disp(tr)