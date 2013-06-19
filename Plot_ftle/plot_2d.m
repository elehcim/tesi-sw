function plot_2d( my_vis_var, my_fixed_var, param_struc, file_content, flags)
n1=param_struc.(my_vis_var{1});
n2=param_struc.(my_vis_var{2});

coord1_grid=reshape(file_content(:,1),n2,n1);
coord2_grid=reshape(file_content(:,2),n2,n1);
ftle=file_content(:,3);
coord1_vec=coord1_grid(1,:);
coord2_vec=coord2_grid(:,1);
ftle_matrix=reshape(ftle,n2,n1);

if flags.save_mat
	v=ftle_matrix;
	X=coord1_vec;
	Y=coord2_vec;
	[~,name,~]=fileparts(param_struc.file_name);
	data_name=sprintf('DATA2d_%s.mat',name);
	uisave({'v','X','Y','param_struc','my_fixed_var','my_vis_var'},...
		[param_struc.file_dir,data_name])
end

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
% TODO insert flags.ftle

%% Remove zeros
ftle_matrix(ftle_matrix==0)=NaN;


%% Create labels
labels = create_labels( my_vis_var );
%% Write generic part of the title
title_string=(sprintf('%ix%i\n%s=%.3f %s=%.3f ecc=%.2f t=%.2f T=%.2f',...
	n1,n2,my_fixed_var{1},param_struc.(my_fixed_var{1}),...
	my_fixed_var{2},param_struc.(my_fixed_var{2}),...
	param_struc.ecc, param_struc.t0, param_struc.DT));
%% Plot data
ftle_fig_handle=figure;
pcolor(coord1_vec, coord2_vec, ftle_matrix)
colorbar; axis square; shading flat; xlabel(labels(1)); ylabel(labels(2));

if isfield(param_struc,'method') && strcmp(param_struc.method,'''FILE''')
	title([sprintf('FILE %i intersection ',param_struc.n_iterations) title_string]);
else
	title(['FTLE ' title_string]);
end
%% Set font size
font_size=20;
set(findall(ftle_fig_handle,'-property','FontSize'),'FontSize',font_size)
%% plot gridfit
if flags.gridfit
	addpath(genpath('gridfitdir'));
	x=file_content(:,1);
	y=file_content(:,2);
	z=file_content(:,3);
	g=gridfit(x,y,z,coord1_vec,coord2_vec);
	figure
	surf(coord1_vec,coord2_vec,g);
	shading interp
end
%% Locate and plot lcs
if flags.lcs
	% old way to compute lcs
	% lcs=locate_lcs(ftle,0.64);
	% lcs_matrix=reshape(lcs,n2,n1);
	% Plot lcs
	lcs_fig_handle=figure;
	pcolor(coord1_vec, coord2_vec, lcs_matrix)
	colorbar; shading flat; xlabel(labels(1)); ylabel(labels(2));
	title(['LCS ' title_string]);
end
%% Save file
if flags.save_fig
		[~,save_fig_name,~]=fileparts(param_struc.file_name);
		set(ftle_fig_handle,'UserData',param_struc)
		hgsave(ftle_fig_handle,[flags.output_dir,save_fig_name,'.fig'])
	if flags.lcs
		error('not yet implemented')
	end
end

