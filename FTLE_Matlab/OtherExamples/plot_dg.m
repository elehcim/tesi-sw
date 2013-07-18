font_size=20;

load data_dg_n300-T5.mat
figure('name','Various double gyre integration time')
%subplot(2,2,1)
pcolor(grid_x,grid_y,sigma');
shading interp
%colorbar('location','EastOutside');
axis image
xlabel('$x$','fontsize',font_size,'interpreter','latex')
ylabel('$y$','fontsize',font_size,'interpreter','latex')
set(gca,'FontSize',font_size)

load data_dg_n300-T10.mat
figure
%subplot(2,2,2)
pcolor(grid_x,grid_y,sigma');
shading interp
%colorbar('location','EastOutside');
axis image
xlabel('$x$','fontsize',font_size,'interpreter','latex')
ylabel('$y$','fontsize',font_size,'interpreter','latex')
set(gca,'FontSize',font_size)

load data_dg_n300-T15.mat
figure
%subplot(2,2,3)
pcolor(grid_x,grid_y,sigma');
shading interp
%colorbar('location','EastOutside');
axis image
xlabel('$x$','fontsize',font_size,'interpreter','latex')
ylabel('$y$','fontsize',font_size,'interpreter','latex')
set(gca,'FontSize',font_size)

load data_dg_n300-T20.mat
figure
%subplot(2,2,4)
pcolor(grid_x,grid_y,sigma');
shading interp
%colorbar('location','EastOutside');
axis image
xlabel('$x$','fontsize',font_size,'interpreter','latex')
ylabel('$y$','fontsize',font_size,'interpreter','latex')
set(gca,'FontSize',font_size)

load data_dg_n300-T-15.mat
figure
%subplot(2,2,4)
pcolor(grid_x,grid_y,sigma');
shading interp
%colorbar('location','EastOutside');
axis image
xlabel('$x$','fontsize',font_size,'interpreter','latex')
ylabel('$y$','fontsize',font_size,'interpreter','latex')
set(gca,'FontSize',font_size)