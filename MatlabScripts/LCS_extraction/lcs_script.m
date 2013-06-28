clearvars
%% 500
%load('data2d_mu=0.1000_ecc=0.00nx=500y0=0.00nvx=500e0=-1.81t0=0.00_1216.mat')
%% Sun Jupiter
%load('data/DATA2d_Sun_Jupiter_4it_bw.mat');
%% 1000
%load('data/data2d_mu=0.1000_ecc=0.00nx=1000y0=0.00nvx=1000e0=-1.81t0=0.00_old.mat')
%% Occhio Gawlik
load('data/DATA2d_OcchioGawlik.mat')
%% DATA
epsilon=1e-4;
resample=0;
filter=1; % 1:gauss; 2:unsharp
sigma=1;
siz=sigma*6;
v(v==0)=NaN;
%% Plot original image
orig=subplot(1,2,1);
pcolor(X,Y,v);shading interp;title('Original Image')
%% Filter
G_gauss = fspecial('gaussian',[siz siz],sigma);
G_unsharp = fspecial('unsharp');
switch filter
	case 1
		v = imfilter(v,G_gauss);
		% v=imgaussian(v,sigma,siz);
	case 2
		v = imfilter(v,G_unsharp);
end

%%
% subplot(2,4,3:4)
% disp('filtered')
% pcolor(X,Y,v);shading interp;colorbar;title('Filtered Image')

ridges=subplot(1,2,2);
LCS2d_plot_test(v,epsilon,X,Y,'cross product',resample);

% subplot(2,4,5:6)
% disp('method: cross product')
% LCS2d_plot_test(v,epsilon,X,Y,'cross product',resample);

% subplot(2,2,3)
% disp('method: curvature vector')
% LCS2d_plot_test(v,epsilon,X,Y,'curvature vector',resample);

subplot(2,4,7:8)
disp('method: eigen vector')
LCS2d_plot_test(v,epsilon,X,Y,'eigen vector',resample);

%% Set font size
h = gcf;
font_size=20;
set(findall(h,'-property','FontSize'),'FontSize',font_size)