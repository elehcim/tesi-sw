function lcs=LCS_plot(Volume, Sigma, epsilon,X,Y,Z,smooth)

%Volume(Volume==0)=nan;

if nargin<7
	smooth =0;
end
if smooth
	[Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(Volume,Sigma);
	if(Sigma>0)
		F=imgaussian(Volume,Sigma);
	else
		F=Volume;
	end
else
	[Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = my_Hessian3D(Volume);
	F=Volume;
end

Dx=gradient3(F,'x');
Dy=gradient3(F,'y');
Dz=gradient3(F,'z');

Hgx=Dxx.*Dx + Dxy.*Dy + Dxz.*Dz;
Hgy=Dxy.*Dx + Dyy.*Dy + Dyz.*Dz;
Hgz=Dxz.*Dx + Dyz.*Dy + Dzz.*Dz;

normHg=sqrt(Hgx.^2+Hgy.^2+Hgz.^2);
[Lambda1,~,~]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
normLambdag=sqrt((Lambda1.*Dx).^2+(Lambda1.*Dy).^2+(Lambda1.*Dz).^2);
n_f=5;
%lcs_cell={1:n_f};
%lcs=nan(size(Volume));

% for i=2:n_f
	dimensions=size(Volume);
	lcs=nan(dimensions);
% 	epsilon=10^-i
% Try to eliminate borders
% 	lcs(1,1,:)=0;
% 	lcs(1,:,1)=0;
% 	lcs(:,1,1)=0;
% 	lcs(dimensions(1),:,:)=0;
% 	lcs(:,dimensions(2),:)=0;
% 	lcs(:,:,dimensions(3))=0;
	lcs(isnan(lcs))=0;
	lcs(abs(normHg-normLambdag)<epsilon)=1;
	% max(lcs(:))
	% min(lcs(:))
% 	lcs_cell{i}=lcs;
% end



[X_grid,Y_grid,Z_grid]=meshgrid(X,Y,Z);
% F=TriScatteredInterp(X,Y,Z)
% z_for_plot=linspace(min(Z),max(Z),length(X));

%for i=2:n_f
%	epsilon=10^-i
	figure
isosurface(X_grid,Y_grid,Z_grid,lcs,.8);
% isonormals(X_grid,Y_grid,Z_grid,lcs,pippo)
% set(pippo,'FaceColor','red','EdgeColor','none')

	
	axis square
pbaspect(gca,dimensions)
%title(sprintf('epsilon=%.1e',epsilon))
%end

% lcs(abs(normHg-normLambdag)<epsilon)=1;
% lcs(isnan(lcs))=0;
% max(lcs(:))
% min(lcs(:))
% patch(isosurface(lcs,1))
% sliceomatic(lcs_cell{2})

%{
x=Points(:,1);
y=Points(:,2);
z=Points(:,3);

dx=1;
dy=1;

x_edge=[floor(min(x)):dx:ceil(max(x))];
y_edge=[floor(min(y)):dy:ceil(max(y))];
[X,Y]=meshgrid(x_edge,y_edge);
Z=griddata(x,y,z,X,Y);
% The following line of code is if you use JE's gridfit:
% Z=gridfit(x,y,z,x_edge,y_edge);

%NOW surf and mesh will work...

surf(X,Y,Z)
%mesh(X,Y,Z)
%}
%{
% Alternative to griddata:
F = TriScatteredInterp(x,y,z);
Z= F(X,Y);
%}
end

function filter_lcs(lcs)
index=size(lcs);
for i=1:index(1)
	for j=1:index(1)
		for k=1:index(1)
		continue
		end
	end
end		
end
 function D = gradient3(F,option)
% This function does the same as the default matlab "gradient" function
% but with one direction at the time, less cpu and less memory usage.
%
% Example:
%
% Fx = gradient3(F,'x');

[k,l,m] = size(F);
D  = zeros(size(F),class(F)); 

switch lower(option)
case 'x'
    % Take forward differences on left and right edges
    D(1,:,:) = (F(2,:,:) - F(1,:,:));
    D(k,:,:) = (F(k,:,:) - F(k-1,:,:));
    % Take centered differences on interior points
    D(2:k-1,:,:) = (F(3:k,:,:)-F(1:k-2,:,:))/2;
case 'y'
    D(:,1,:) = (F(:,2,:) - F(:,1,:));
    D(:,l,:) = (F(:,l,:) - F(:,l-1,:));
    D(:,2:l-1,:) = (F(:,3:l,:)-F(:,1:l-2,:))/2;
case 'z'
    D(:,:,1) = (F(:,:,2) - F(:,:,1));
    D(:,:,m) = (F(:,:,m) - F(:,:,m-1));
    D(:,:,2:m-1) = (F(:,:,3:m)-F(:,:,1:m-2))/2;
otherwise
    disp('Unknown option')
end
end
