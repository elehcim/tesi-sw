function [my_vis_var,my_fixed_var,param_struc,file_content]=parse_file(file_name)
fid=fopen(file_name,'r');
if (fid < 0)
	error('Error:could not open file %s',file_name)
else
	line = fgetl(fid);
	nl=1;
	while ~(strcmp(line,'### Beginning of FTLE data ###') ||... % for retrocompatibility
			strcmp(line,'### Beginning of data ###') ||...
			feof(fid))
		tmp=textscan(line,'%s %f','Delimiter','=');
		if isempty(tmp{1,2})
			param{nl}=textscan(line,'%s %s','Delimiter','=');
			param_struc.(param{nl}{1}{1})=param{nl}{2}{1};
		else
			param{nl}=textscan(line,'%s %f','Delimiter','=');
			param_struc.(param{nl}{1}{1})=param{nl}{2};
		end
		
		line=fgetl(fid);
		nl=nl+1;
	end;
end;

%% add other fields to the Parameter structure
% put the file_name in the struct
[file_dir,name,ext]=fileparts(file_name);
param_struc.file_dir=fullfile(file_dir,filesep);
param_struc.file_name=[name,ext];


%% Parse variables to plot and fixed variables
vis_var={'nx','ny','nvx','nvy','ne'};
fixed_var={'x_0','y_0','vx_0','vy_0','e_0'};
c1=1;
c2=1;
my_vis_var=cell(2,1);my_fixed_var=cell(2,1);
for i=1:length(param)
	if any(strcmp(param{i}{1},vis_var))
		my_vis_var{c1}=param{i}{1}{1};
		c1=c1+1;
	elseif any(strcmp(param{i}{1},fixed_var))
		my_fixed_var{c2}=param{i}{1}{1};
		c2=c2+1;
	end
end
%% Select dimension
dim=length(my_vis_var);
if dim==2
	formatSpec='%f\t%f\t%f';
elseif dim==3
	formatSpec='%f\t%f\t%f\t%f';
else
	fclose(fid);
	error('I can''t determine dimensions of the ftle data')
end
%% Read data
ftle_data=textscan(fid,formatSpec);
file_content=cell2mat(ftle_data);
fclose(fid);
end