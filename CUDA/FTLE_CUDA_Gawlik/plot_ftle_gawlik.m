n=150;
filename='ftle_ell_n=150_t=2.09.txt';
fid=fopen(filename,'r');
% %% Determina n
% [~,n]=fscanf(fid,'\n')
% n
%TODO provare questa
%[stat, n_l_s] = system(['grep -c ".$" ' filename]);
%n_lines = str2double(n_l_s)
% o anche questa
% file = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
% n_lines=lenght(file)

%% Initialize axis
nx=n;
x_0_min=-0.8;
x_0_max=-0.15;
dx=(x_0_max-x_0_min)/(nx-1);
x_0=x_0_min : dx : x_0_max;

nvx=n;
vx_0_min=-2;
vx_0_max=2;
dvx=(vx_0_max-vx_0_min)/(nvx-1);
vx_0=vx_0_min : dvx : vx_0_max;

%% Read file
ftle=fscanf(fid,'%g',[n n]);
% size(ftle)
%% Plot
figure
pcolor(x_0, vx_0, ftle)
colorbar
shading flat
fclose(fid);