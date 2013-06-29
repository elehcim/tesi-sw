clear all
% circolare
%file_name='txt/tracers_mu=0.1000_ecc=0.00nx=100y0=0.00nvx=100e0=-1.87.txt';
%file_name='txt/tracers_mu=0.1000_ecc=0.00nx=10y0=0.00nvx=10e0=-1.87.txt';
%file_name='txt/tracers_mu=0.1000_ecc=0.00nx=50y0=0.00nvx=50e0=-1.87.txt';
%file_name='txt/tracers_mu=0.1000_ecc=0.00nx=50y0=0.00nvx=50e0=-1.74.txt';

% ellittico
%file_name='tracers_mu=0.1000_ecc=0.04nx=50y0=0.00nvx=50e0=-1.87.txt';
%file_name='tracers_mu=0.1000_ecc=0.04nx=50y0=0.00nvx=50e0=-1.70.txt';
%file_name='tracers_mu=0.1000_ecc=0.04nx=50y0=0.00nvx=50e0=-1.74.txt';

% file_name='tracers_mu=0.1000_ecc=0.04nx=25y0=0.00nvx=25e0=-1.74_SACRO.txt';

setappdata(0,'UseNativeSystemDialogs',false)
file_name=uigetfile('.txt','','txt/');
file_name=['txt/',file_name;]
[~,name,~]=fileparts(file_name);
mat_name=['mat/', name,'.mat']
if exist(['mat/',mat_name],'file')
	disp('opening .mat file')
	load(mat_name)
else
	disp('parsing .txt file')
	formatSpec='%f\t%f\t%f\t%f\t%f';
	fid=fopen(file_name,'r');
	if (fid < 0)
		error('Error:could not open file %s', file_name)
	else
		time_counter=1;
		while ~feof(fid)
			line = fgets(fid);
			if strcmp(line(1),'#')
				first_cut=strrep(line,'### t=','');
				time_step(time_counter)=str2double(strrep(first_cut,' ###',''));
				time_counter=time_counter+1;
				i=1;
			else
				tracers{time_counter-1,1}(i,:)=cell2mat(textscan(line,formatSpec));
				i=i+1;
			end
		end
	end

n_tracers=i-1;
n_frames=time_counter-1;

% take data from file name
mu=str2double(file_name(12:17))
ecc=str2double(file_name(23:26))
% take data from the file
e=tracers{1,1}(1,end)
save(mat_name)

end
%% Identify which tracers are inside a LCS
% Setting interval of the LCS
% x_lcs_min
% x_lcs_max
% vx_lcs_min
% vx_lcs_max

fprintf('I have %i tracers\n',n_tracers)%  of which %i where inside the LCS\n',n_tracers,4)%FIXME
%% make figure
pause_time=.5; %s

figure('NumberTitle','off');%,'name', sprintf('Interval = %.1fs',pause_time))

%% plot pseudo-Hill's region
points=500;

bb=3; % Bounding box
x=linspace(-bb,bb,points);
y=linspace(-bb,bb,points);
[x,y]=meshgrid(x,y);
%z=(Omega(x,y,mu));
%{
 --- Useless ---
[Hill,h]=contour(x,y,z,[-e,-e]);
x1_hill=Hill(1,2:Hill(2,1)+1);
y1_hill=Hill(2,2:Hill(2,1)+1);
x2_hill=Hill(1,Hill(2,1)+3:end);
y2_hill=Hill(2,Hill(2,1)+3:end);
%}
%% Plot

% select colour
% load('lcs_points_1')
load inner_spline
xys_inner=xys;
load outer_spline
xys_outer=xys;

ecc=0.04

%ecc=0.04839;
%mu=9.537e-4 Sun Jup
mu=0.1;
c=nan(1,n_tracers);
d=nan(1,n_tracers);
for k=1:n_tracers
	% I want to plot the inner tracers
	% I want also to remove the tracers in between the inner and the outer splines
	c(k)=inpolygon(tracers{1}(k,1),tracers{1}(k,3),xys_inner(1,:),xys_inner(2,:));
	d(k)=inpolygon(tracers{1}(k,1),tracers{1}(k,3),xys_outer(1,:),xys_outer(2,:));
end
dd=d & ~c;
for k=1:n_frames
	tracers{k}(dd,:)=[];
end
c(dd)=[];
size(tracers{1})
size(c)
% plot
for j=[1 7 20 31]%j=1:n_frames
 	figure
	x_t=tracers{j}(:,1);
	y_t=tracers{j}(:,2);
	scatter(x_t,y_t,20,c,'o','filled')
	hold on
 	z=Omega(x,y,mu)/(1+ecc*cos(time_step(j)));
 	contour(x,y,z,[-e,-e]);
	xlim([-1.2 1.2])
	ylim([-1.2 1.2])
	xlabel('$x$','fontsize',13,'interpreter','latex')
	ylabel('$y$','fontsize',13,'interpreter','latex')
	hold off
	%title(sprintf('ecc=%.2f, e_0=%.2f, t=%.2f, frame %i of %i',ecc,e,time_step(j),j,n_frames))
	title(sprintf('t=%.2f',time_step(j)))
	pause(pause_time);
end
