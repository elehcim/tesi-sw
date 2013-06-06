function plot_3d( my_vis_var, my_fixed_var, param_struc, file_content, flags)
n1=param_struc.(my_vis_var{1});
n2=param_struc.(my_vis_var{2});
n3=param_struc.(my_vis_var{3});
% Check if the file is coherent

% adjust
coord1_grid=reshape(file_content(:,1),n3,n2,n1);
coord2_grid=reshape(file_content(:,2),n3,n2,n1);
coord3_grid=reshape(file_content(:,3),n3,n2,n1);
coord1_vec=squeeze(coord1_grid(1,1,:));
coord2_vec=squeeze(coord2_grid(1,:,1))';
coord3_vec=squeeze(coord3_grid(:,1,1));
ftle=file_content(:,4);
ftle_matrix=permute(reshape(ftle,n3,n2,n1),[2 3 1]);
% ftle_matrix=squeeze(reshape(ftle,n3,n2,n1));

%% Create labels
labels = create_labels( my_vis_var );
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

%% Save mat
if flags.save_mat
	v=ftle_matrix;
	X=coord1_vec;
	Y=coord2_vec;
	Z=coord3_vec;
	[~,name,~]=fileparts(param_struc.file_name);
	data_name=sprintf('DATA3d_%s.mat',name);
	%data_name=sprintf('data3d_%s.mat',param_struc.file_name(10:end-4));
	uisave({'v','X','Y','Z','param_struc','my_fixed_var','my_vis_var'},...
		[param_struc.file_dir,data_name])
end

%% Locate lcs %FIXME
% 
% lcs_matrix=locate_lcs(ftle_matrix,0.6);

%% Write generic part of the title
title_string=(sprintf('%ix%ix%i\n%s=%.3f ecc=%.2f t=%.2f T=%.2f',...
	n1,n2,n3,my_fixed_var{1},param_struc.(my_fixed_var{1}),...
	param_struc.ecc, param_struc.t0, param_struc.DT));
%FIXME
%% Plot ftle
% f=figure;
% 	[x1,x2,x3]=meshgrid(coord1_vec, coord2_vec, coord3_vec);
% 	slice(x1,x2,x3,ftle_matrix,[],[],-1.7);
% colorbar; shading flat; xlabel(labels(1)); ylabel(labels(2)); zlabel(labels(3));
% title(['FTLE ' title_string]);xlabel(labels(1)); ylabel(labels(2)); zlabel(labels(3));

%% Plot lcs
% g=figure;
% 	[x1,x2,x3]=meshgrid(coord1_vec, coord2_vec, coord3_vec);
% 	slice(x1,x2,x3,lcs_matrix,[],[],-1.7);
% colorbar; shading flat; xlabel(labels(1)); ylabel(labels(2)); zlabel(labels(3));
% title(['LCS ' title_string]);xlabel(labels(1)); ylabel(labels(2)); zlabel(labels(3));

%% myslicer
% fprintf('\n\nMove the slices with your mouse pointer. If you use any of the standard tools that \ntakes control of the mouse (e.g. rotation), be sure to click on them again in order to \ndeactivate them before you try to move the slices again.\n\n')
% 
% figure;
% camproj('perspective')
% T = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
% myslicer(squeeze(ftle_matrix),T);
% %colormap gray;
% %axis off;
% colorbar; xlabel(labels(1)); ylabel(labels(2)); zlabel(labels(3));

%% Sliceomatic FTLE
% sliceomatic(ftle_matrix)
addpath(genpath([folder 'Plot_library/sliceomatic-2.3']));
sliceomatic(ftle_matrix,coord1_vec,coord2_vec,coord3_vec)
axis square
pbaspect(gca,[n1 n2 n3])
title(['FTLE ' title_string]);
xlabel(labels(1));ylabel(labels(2));zlabel(labels(3));
%% Set font size
h = gcf;
font_size=12;
set(findall(h,'-property','FontSize'),'FontSize',font_size)
%% Sliceomatic LCS
% sliceomatic(lcs_matrix)
% title(['LCS ' title_string]);xlabel(labels(1)); ylabel(labels(2)); zlabel(labels(3));

end

