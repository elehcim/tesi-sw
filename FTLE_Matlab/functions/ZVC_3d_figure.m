%% Method 1: ezimplot3
f=@(x,y,z)(x^2)/2+(y^2)/2+0.9/(sqrt((x+0.012)^2+y^2+z^2))+0.1/(sqrt(x-1+0.012)^2+y^2+z^2)-1.45;
figure
ezimplot3(f,[-1.5,1.5,-1.5,1.5,-1.5,1.5])

%% Method 2: isosurface
[x,y,z]=meshgrid(-1.5:0.01:1.5);
F=(x.^2)/2+(y.^2)/2+0.9./(sqrt((x+0.012).^2+y.^2+z.^2))+0.1./(sqrt(x-1+0.012).^2+y.^2+z.^2);
figure
isosurface(x,y,z,F,1.70)