function plot_tracers(tracers,contours)
%% select colour
length(contours)

% load('lcs_points_1')
load inner_spline
xys_inner=xys;
load outer_spline
xys_outer=xys;

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
for j=1:n_frames
	x_t=tracers{j}(:,1);
	y_t=tracers{j}(:,2);
	scatter(x_t,y_t,20,c,'o','filled')
	hold on
	z=Omega(x,y,mu)/(1+ecc*cos(time_step(j)));
	contour(x,y,z,[-e,-e]);
	hold off
	title(sprintf('ecc=%.2f, e_0=%.2f, t=%.2f, frame %i of %i',ecc,e,time_step(j),j,n_frames))
	pause(pause_time);
end
