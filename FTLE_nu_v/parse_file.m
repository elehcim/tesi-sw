function [param_struc,file_content]=parse_file(file_name)
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
		if isempty(tmp{1,2}) % If it's a string
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

formatSpec='%f\t%f\t%f';

%% Read data
ftle_data=textscan(fid,formatSpec);
file_content=cell2mat(ftle_data);
fclose(fid);
end