clear all
%% Locations of folders and files
gif_name='Ellip200.gif';
gif_dir=[folder,'ER3BP/Figure/Video/',gif_name];
pic_dir=[folder,'ER3BP/Figure/Figure.6/'];
n_frames=20;

%% Preallocation
pic_list=cell(n_frames,1);
% Preallocate movie structure.
% (copied form "movie2avi" help)
M(1:n_frames) = struct('cdata', [], 'colormap', []);

%% Create a list of figure names
listing=dir(pic_dir);
for i=1:n_frames
	pic_list{i}=[pic_dir listing(i+2).name];
end

%% Load figures in an array of handles
for j=1:n_frames
	%h(j)=openfig(pic_list{j},'reuse','invisible');
	h(j)=openfig(pic_list{j});
end
%set(h(:),'Visible','off','Renderer','zbuffer','NextPlot','replacechildren');
set(h(:),'Position',[100 100 850 600],'Renderer','zbuffer') % TODO zbuffer seems the best

%% Build the movie for FTLE

for j=1:n_frames % era 1:10
	M(j)=getframe(h(j)); % era h(j+10)
	% pause(0.5); % seems needed to stabilize the frames
	% ---- start test for GIF ----
	% TODO set the interval
	% adapted from: http://www.mathworks.it/support/solutions/en/data/1-48KECO/
	im = frame2im(M(j));
 	[imind,cm] = rgb2ind(im,256);
	if j == 1;
		imwrite(imind,cm,gif_dir,'gif', 'Loopcount',inf);
	else
		imwrite(imind,cm,gif_dir','gif','WriteMode','append');
	end
	% ---- end test for GIF ----
end
close all % close all the getframe figures
%% Adjust properties of the frames
% Copied from "movie" help
% use 1st frame to get dimensions
[height, w, p] = size(M(1).cdata);
hf=figure;
% resize figure based on frame's w x h, and place at (150, 150)
set(hf,'Position', [150 150 w height]);
axis off
%% Display movie
times=-5;
fps=12;
movie(hf,M,times,fps)
%% tests for movie
% movie2avi(M, '../Video/Ellip.avi', 'compression', 'None', 'quality', 100);
% TODO next step oto obtain avi is to try VideoWriter