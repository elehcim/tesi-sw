function plot_multiple_files(data_dir,output_dir)
% Plots files in data_dir and save it in out_dir
listing=dir(data_dir);
data_list=cell(size(listing,1)-2,1);
disp('Data files:')
number_of_files=size(listing,1)-2;
for i=1:number_of_files
	data_list{i}=[data_dir listing(i+2).name];
	disp(data_list{i}(length(data_dir)+1:end));
end
%% Read n_frames from first file
n_frames=parse_n_frames(data_list{1});
if n_frames ~= number_of_files
	fprintf('n_frames ~= number of files.\nTaking present files\n')
	n_frames=number_of_files;
end
%% Read data and write figures
for i=1:n_frames
	%% Read file
	current_file=data_list{i}(length(data_dir)+1:end);
	fprintf('Analyzing file %i of %i: \t%s\n',i,n_frames,current_file);
	[my_vis_var,my_fixed_var,param_struc,file_content]=...
		parse_file([data_dir current_file]);
	n1=param_struc.(my_vis_var{1});
	n2=param_struc.(my_vis_var{2});
	
	coord1_grid=reshape(file_content(:,1),n2,n1);
	coord2_grid=reshape(file_content(:,2),n2,n1);
	ftle=file_content(:,3);
	coord1_vec=coord1_grid(1,:);
	coord2_vec=coord2_grid(:,1);
	ftle_matrix=reshape(ftle,n2,n1);
	
	% locate lcs
% 	lcs=locate_lcs(ftle,0.64);
% 	lcs_matrix=reshape(lcs,n2,n1);
	%% Create labels
	labels = create_labels( my_vis_var );
	% plot ftle
	ftle_fig=figure;
	pcolor(coord1_vec, coord2_vec, ftle_matrix)
	colorbar; shading flat; xlabel(labels(1)); ylabel(labels(2));
	title(sprintf('FTLE %ix%i\n%s=%.3f %s=%.3f ecc=%.2f t=%.2f',...
		n1,n2,my_fixed_var{1},param_struc.(my_fixed_var{1}),...
		my_fixed_var{2},param_struc.(my_fixed_var{2}),...
		param_struc.ecc, param_struc.t0));
	
	% plot lcs
% 	lcs_fig=figure;
% 	pcolor(coord1_vec, coord2_vec, lcs_matrix)
% 	colorbar; shading flat; xlabel(labels(1)); ylabel(labels(2));
% 	title(sprintf('LCS %ix%i\n%s=%.3f %s=%.3f ecc=%.2f t=%.2f',...
% 		n1,n2,my_fixed_var{1}{1},param_struc.(my_fixed_var{1}{1}),...
% 		my_fixed_var{2}{1},param_struc.(my_fixed_var{2}{1}),...
% 		param_struc.ecc, param_struc.t0));
	% write
	save_plots(ftle_fig,output_dir,n1,n2,param_struc.t0)
end