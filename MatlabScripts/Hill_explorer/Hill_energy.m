clear all
 points=500;

mu=9.537e-4;
e_min=-2;
e_max=-1;
ne=100;
bb=2; % Bounding box
x=linspace(-bb,bb,points);
y=linspace(-bb,bb,points);
e=linspace(e_min,e_max,ne);
[xx,yy]=meshgrid(x,y);
EE=Omega(xx,yy,mu)


%TODO

% clear all
% %% plot pseudo-Hill's region
% points=500;
% n_frames=100;
% mu=9.537e-4;
% e_min=-2;
% e_max=-1;
% energy_vector=linspace(e_min,e_max,n_frames);
% pause_time=.1;
% bb=2; % Bounding box
% x=linspace(-bb,bb,points);
% y=linspace(-bb,bb,points);
% [x,y]=meshgrid(x,y);
% j=1;
% figure
% for e=energy_vector
% 	z=Omega(x,y,mu);
% 	contour(x,y,z,[-e,-e]);
% 	title(sprintf('e_0=%.2f, frame %i of %i',e,j,n_frames))
% 	pause(pause_time);
% 	j=j+1;
% end