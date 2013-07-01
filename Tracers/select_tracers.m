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
			('*.fig','Select FTLE figure');
	end
	fig_file=fullfile(PathName,fig_file);
end

[~,fig_file_name,~]=fileparts(fig_file);
fig=openfig(fig_file);

%% Get the param_struc structure in the figure UserData field.
set(gcf, 'Renderer', 'zbuffer');
param_struc=get(fig,'UserData');
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
	plot(xi,yi,'+k','MarkerSize',20)
	n = n+1;
	fprintf('tracer %i: %.2f\t%.2f\n',n,xi,yi)
	text(xi,yi,sprintf(' %d',n),'fontsize',25,'VerticalAlignment','top')
	xy(:,n) = [xi;yi];
end
hold off

%% Add missing coordinates
vars={'x','y','vx','vy','e'};
nvars={'nx','ny','nvx','nvy','ne'};
vars_0={'x_0','y_0','vx_0','vy_0','e_0'};
tr.n_tracers=n;
tr.x=nan(1,n);
tr.y=nan(1,n);
tr.vx=nan(1,n);
tr.vy=nan(1,n);
tr.e=nan(1,n);
c=1;
for i=1:5
	if isfield(param_struc,nvars{i})
		tr.(vars{i})=xy(c,:);
		c=c+1;
	elseif isfield(param_struc,vars_0{i})
		tr.(vars{i})=param_struc.(vars_0{i})*ones(1,n);
	end
end

%% Add other parameters
tr.t0=param_struc.t0;
tr.T=param_struc.DT;
tr.mu=param_struc.mu;
tr.ecc=param_struc.ecc;

%% Reconstruct missing coordinate
tr=complete_tracers(tr);

%% save data
time_stamp=clock;
clo=sprintf('%04d%02d%02d-%02.0f%02.0f%02.0f',time_stamp(1),time_stamp(2),...
	time_stamp(3),time_stamp(4),time_stamp(5),time_stamp(6));
data_name=[fig_file_name '_tracers_' clo '.mat'];
save(data_name,'tr')