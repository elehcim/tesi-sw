function [tr]=select_tracers(varargin)
% Select tracers from a fig file.
% It's possible to pass a string containing the path to a fig file that 
% will be the input.
% Can accept a directory that will be the working directory for uigetfile.
%% Parse input
if length(varargin)>1
	error('too many arguments')
elseif length(varargin)==1
	if exist(varargin{1},'dir');
		working_dir=varargin{1};
	elseif exist(varargin{1},'file')
		fig_file=varargin{1};
	end
end
%% Open figure if it wasn't passed to the function
if ~exist('fig_file','var')
	if length(varargin)==1
		fprintf('Can''t find file: %s\nPlease specific one .fig file', varargin{1});
	end
	if exist('working_dir','var')
		[fig_file, PathName]=uigetfile...
			('*.fig','Select FTLE figure',working_dir);
	else
		[fig_file, PathName]=uigetfile...
			('*.fig','Select FTLE figure',[pwd '/../']);
	end
	fig_file=fullfile(PathName,fig_file);
end

[~,fig_file_name,~]=fileparts(fig_file);
fig=openfig(fig_file);

%% Get the param_struc structure in the figure UserData field.
set(gcf, 'Renderer', 'zbuffer');
param_struc=get(fig,'UserData')
% %% get image coordinates
% h=get(gca);
% coord_x=get(h.XLabel,'String');
% coord_y=get(h.YLabel,'String');
%% Select tracers
hold on
% Initially, the list of points is empty.
xy = [];
n = 0;
% Loop, picking up the points.
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1;
while but == 1
	[xi,yi,but] = ginput(1);
	plot(xi,yi,'.w','MarkerSize',25)
	n = n+1;
	fprintf('tracer %i: %.2f\t%.2f\n',n,xi,yi)
	text(xi,yi,sprintf(' %d',n),'color',[1 1 1],'fontsize',25,'FontWeight','bold','VerticalAlignment','top')
	% bold 25 white pallino 25
	xy(:,n) = [xi;yi];
end
hold off

%% Add other parameters
tr.n_tracers=n;
tr.t0=param_struc.t0;
tr.T=param_struc.DT;
tr.mu=param_struc.mu;
tr.ecc=param_struc.ecc;

%% Compute earth velocity
a_jup=778412027; %km
GM_jup=126711995; %km^3/s^2
GM_sun=132712439935; %km^3/s^2
GM=GM_jup+GM_sun;
v_E=29.783; % km/s
v_tilde=sqrt(GM*(1+tr.ecc*cos(tr.t0))/(a_jup*(1-tr.ecc^2)));
v_E_adim=v_E/v_tilde;

%% compute initial points
dist_Sun_Jup = 778547200;	% km
r_earth_orbit = 149600000;	% km
R = r_earth_orbit/dist_Sun_Jup;
nu=xy(1,:);

% Add Earth velocity to v
v=xy(2,:)+v_E_adim;

[tr.x,tr.y]=nu2xy(nu,R);
for i=1:n
	tr.vx(1,i)=-v(i)*sin(nu(i));
	tr.vy(1,i)=v(i)*cos(nu(i));
end


%% save data
time_stamp=clock;
clo=sprintf('%04d%02d%02d-%02.0f%02.0f%02.0f',time_stamp(1),time_stamp(2),...
	time_stamp(3),time_stamp(4),time_stamp(5),time_stamp(6));
data_name=[fig_file_name '_tracers_' clo '.mat'];
save(data_name,'tr')