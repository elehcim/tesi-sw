function plot_multiple_files(files, data_dir, flags)
number_of_files=length(files);
%% List file
disp('Data files:')
for i=1:number_of_files
	disp(files{i});
end
%% Read n_frames from first file
%if it's not coherent with the number of file use the number of files
path_to_first_file=fullfile(data_dir,files{1});
n_frames=parse_n_frames(path_to_first_file);
if n_frames ~= number_of_files
	fprintf('n_frames ~= number of files.\nTaking present files\n')
	n_frames=number_of_files;
end
%% Read data and write figures
h=zeros(n_frames,1); % preallocate
for i=1:n_frames
	%% Read file
	current_file=files{i};
	fprintf('Analyzing file %i of %i: \t%s\n',i,n_frames,current_file);
	[my_vis_var,my_fixed_var,param_struc,file_content]=...
		parse_file(fullfile(data_dir,current_file));
	n1=param_struc.(my_vis_var{1});
	n2=param_struc.(my_vis_var{2});
	
	coord1_grid=reshape(file_content(:,1),n2,n1);
	coord2_grid=reshape(file_content(:,2),n2,n1);
	ftle=file_content(:,3);
	coord1_vec=coord1_grid(1,:);
	coord2_vec=coord2_grid(:,1);
	ftle_matrix=reshape(ftle,n2,n1);
    
    %% Remove spikes
    if flags.spikes
        %check if it's useless
        if max(ftle_matrix(:))<flags.spikes_value
            warning('no spikes removed')
        else
            %TODO better way to identify out of scale points
            ftle_matrix(ftle_matrix>flags.spikes_value)=0;
        end
    end
    
	% locate lcs
	% 	lcs=locate_lcs(ftle,0.64);
	% 	lcs_matrix=reshape(lcs,n2,n1);
	%% Create labels
	labels = create_labels( my_vis_var );
	% plot ftle
	h(i)=figure('Color',[1 1 1]);
	pcolor(coord1_vec, coord2_vec, ftle_matrix)
	colorbar; shading flat; xlabel(labels(1)); ylabel(labels(2));
	title(sprintf('FTLE %ix%i\n%s=%.3f %s=%.3f ecc=%.2f t=%.2f T=%.2f',...
		n1,n2,my_fixed_var{1},param_struc.(my_fixed_var{1}),...
		my_fixed_var{2},param_struc.(my_fixed_var{2}),...
		param_struc.ecc, param_struc.t0, param_struc.DT));
	%% Set font size
	font_size=20;
	set(findall(h(i),'-property','FontSize'),'FontSize',font_size)
	% % plot lcs
	% 	lcs_fig_handle=figure;
	% 	pcolor(coord1_vec, coord2_vec, lcs_matrix)
	% 	colorbar; shading flat; xlabel(labels(1)); ylabel(labels(2));
	% 	title(sprintf('LCS %ix%i\n%s=%.3f %s=%.3f ecc=%.2f t=%.2f',...
	% 		n1,n2,my_fixed_var{1}{1},param_struc.(my_fixed_var{1}{1}),...
	% 		my_fixed_var{2}{1},param_struc.(my_fixed_var{2}{1}),...
	% 		param_struc.ecc, param_struc.t0));
	
	% write
	if flags.save_fig
		% remove extension and save
		[~,save_fig_name,~]=fileparts(param_struc.file_name);
		hgsave(h(i),[flags.output_dir,save_fig_name,'.fig'])
	end
end
%% Make gif
if flags.gif
	%% build the file name
	
	[~,fname,~]=fileparts(files{1});
	gif_name=[fname '.gif'];
	%% set size
	gif_size=[700 600]; % pixels, please take into account the colorbar
	set(h(:),'Position',[100 100 gif_size],'Renderer','zbuffer');
	
	%% Build the movie for FTLE
	M(1:n_frames) = struct('cdata', [], 'colormap', []); % preallocate
	for j=1:n_frames
		M(j)=getframe(h(j));
		% TODO set the interval
		% adapted from: http://www.mathworks.it/support/solutions/en/data/1-48KECO/
		im = frame2im(M(j));
		[imind,cm] = rgb2ind(im,256);
		if j == 1;
			imwrite(imind,cm,[flags.output_dir,gif_name],'gif', 'Loopcount',inf);
		else
			imwrite(imind,cm,[flags.output_dir,gif_name]','gif','WriteMode','append');
		end
	end
end
