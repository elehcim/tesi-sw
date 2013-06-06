clear all
%% plot pseudo-Hill's region
points=500;
n_frames=50;
mu=9.537e-4;
ecc=0.048775;
t_0=0;
t_f=0;

handle = uicontrol('Name',Value,...)
if t_0 == t_f, n_frames=1; end
%{
e=-1.53
pause_time=.1;
bb=2; % Bounding box
x=linspace(-bb,bb,points);
y=linspace(-bb,bb,points);
[x,y]=meshgrid(x,y);

figure
surfc(x,y,-log(omega),'EdgeColor','none')
figure
contour(x,y,-log(omega),-0.611)
axis equal

% z=Omega(x,y,mu);
%{
 --- Useless ---
[Hill,h]=contour(x,y,z,[-e,-e]);
x1_hill=Hill(1,2:Hill(2,1)+1);
y1_hill=Hill(2,2:Hill(2,1)+1);
x2_hill=Hill(1,Hill(2,1)+3:end);
y2_hill=Hill(2,Hill(2,1)+3:end);
%}
time_step=linspace(t_0,t_f,n_frames);
figure
xlim('auto')
for j=1:n_frames
	z=Omega(x,y,mu)/(1+ecc*cos(time_step(j)));
	contour(x,y,z,[-e,-e]);
	title(sprintf('ecc=%.2f, e_0=%.2f, t=%.2f, frame %i of %i',ecc,e,time_step(j),j,n_frames))
	pause(pause_time);
end
%}