function [contours]=select_contours(varargin)
% select contours
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
end
[~,fig_file_name,~]=fileparts(fig_file);
fig=openfig(fig_file);


%TODO Allow for multiple contours
% prompt={'Multi contours case: number of contours'};
% dlg_title='number of contours';
% defaultanswer={'5'};
% answer=cell2mat(inputdlg(prompt,dlg_title,1,defaultanswer))

%% Select contours
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
	plot(xi,yi,'bo')
	n = n+1;
	xy(:,n) = [xi;yi];
end
% Interpolate with a spline curve and finer spacing.
finer_spacing=.1;
t = 1:n;
ts = 1: finer_spacing : n;
xys = spline(t,xy,ts);

% Plot the interpolated curve.
plot(xys(1,:),xys(2,:),'b-');
hold off
contours=xys;

%% save data
time_stamp=clock;
clo=sprintf('%.0f:%.0f:%.0f',time_stamp(4),time_stamp(5),time_stamp(6));
data_name=[fig_file_name '_contours_' clo];
save(data_name,'xy','xys')

