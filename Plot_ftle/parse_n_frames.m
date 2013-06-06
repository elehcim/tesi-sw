function n_frames=parse_n_frames(path_to_file)
fid=fopen(path_to_file,'r');
if (fid < 0)
	[~,file_name,~]=fileparts(path_to_file);
	error('Error:could not open file %s', file_name)
else
	line = fgetl(fid);
	nl=1;
	param{2}{1}=0;
	while ~feof(fid)
		param{nl}=textscan(line,'%s %f','Delimiter','=');
		line=fgetl(fid);
	if strcmp(param{nl}{1},'n_frames')
		n_frames=param{nl}{2};
		return
	end
	nl=nl+1;
	end
	error('I can''t find n_frames parameter')
end;