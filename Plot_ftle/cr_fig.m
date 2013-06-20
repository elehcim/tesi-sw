function cr_fig ( varargin )
% cr_fig plots FTLE and LCS on .fig files in the provided output dir.
% Input parameters are:
%
% * 'file'
% * 'dir'
% * 'out_dir'
% * 'working_dir'
% 
% 'out_dir' is the output directory where things will be saved.
% If there's only 'file', plot that file in 'out_dir'.
% If there's only 'dir', plot and save files in that directory in 'out_dir'.
% If both are present, plot that single file which is in that directory
% 'dir' in 'out_dir'.
% 'working_dir' is an optional parameter to set the workspace folder, where
% there are data.

%% Count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
	error('cr_fig needs (propertyName , propertyValue) pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
	if any(strcmpi(pair{1},'file'))
		file_name = pair{2};
	elseif any(strcmpi(pair{1},'dir'))
		data_dir = pair{2};
	elseif any(strcmpi(pair{1},'out_dir'))
		output_dir = pair{2};
	elseif any(strcmpi(pair{1},'working_dir'))
		working_dir = pair{2};
	end
end

%% Set default working_dir
if ~exist('working_dir','var') % 
	% FIXME 
	% if on other computers throws an error about 'folder' fulction
	% add this condition to the if clause"if ... && exist('folder','file')"
	
	% TODO verify if directory names end up with a '/'. If not raise an
	% exception
	
	% Set user defaults
	if ispc
		%working_dir=[folder '../FTLE_Calculator/Examples/'];
		working_dir=[folder 'FTLE_CPU/bin/WIN/'];
	else
		% working_dir=[folder 'C++/FTLE_CPU/bin/Linux/'];
		% working_dir=[folder 'FTLE_CPU/bin/Linux/'];
		working_dir=[folder 'MissionAnalysis/Prove per missione/9luglio/']; 
		%working_dir=['/home/michele/Scrivania/1-giu/'];
		% working_dir=[folder 'FTLE_CPU/bin/WIN/'];
		% working_dir='/home/michele/Scrivania/tests/uitest/';
	end
	
	% If no defaults have been yet set, use the current directory
	if ~exist('working_dir','var'), working_dir=pwd; end
end

%% Load file in case it's not provided as input
if ~exist('file_name','var') %usejava('Desktop')
	[FileName, PathName]=uigetfile...
		('*.txt','Select FTLE file',working_dir,'MultiSelect', 'on');
	if length(FileName)==1 && FileName == 0
		return
	else
		file_name=FileName;
		data_dir=PathName;
	end
else
	data_dir=pwd;
end
%% Set flags
flags=opts;
%TODO Set flags also by input to cr_fig()
% example: cr_fig('save_mat',1) should save mat
%%
if ~exist('output_dir','var')
	if flags.save_fig || flags.gif
		out_dir=uigetdir(working_dir,'Select a folder where to save things');
		if out_dir==0
			return
		end
		flags.output_dir=[out_dir, filesep];
	end
end
%% Parse file
if ~iscell(file_name)
	[~,name,ext]=fileparts(file_name);
	fprintf('Loading file: %s\n',[name,ext]);
	path_to_file=fullfile(data_dir,file_name);
	[my_vis_var,my_fixed_var,param_struc,file_content]=parse_file(path_to_file);
	%% Plot single file
	switch length(my_vis_var)
		case 2
			plot_2d( my_vis_var, my_fixed_var, param_struc, ...
				file_content, flags)
		case 3
			plot_3d( my_vis_var, my_fixed_var, param_struc, ...
				file_content, flags)
	end
	%% Plot multiple files
else
	plot_multiple_files(file_name, data_dir, flags)
end
end