function lcs=LCS2d_plot(Image, epsilon,X,Y)
%function lcs=LCS2d_plot(Image, Sigma, ker, epsilon,X,Y)
% if nargin < 2, Sigma = 1;
% elseif nargin < 3 , epsilon = 0.001;
% end
%
% if(Sigma>0)
%     Image=imgaussian(Image,Sigma,ker);
% end
%% Compute Gradient and Hessian
grad_comp='mne';
switch grad_comp
	case 'mine'
		Dx=gradient2(Image,'x');
		Dy=gradient2(Image,'y');
	otherwise
		[Dx,Dy]=gradient(Image);
end
% %% Try to plot gradient
% [X_grid,Y_grid]=meshgrid(X,Y,[0 1]);
% coneplot(X_grid,Y_grid,zeros(size(X_grid)),Dx,Dy,zeros(size(X_grid)),X_grid,Y_grid,zeros(size(X_grid)))

%% Hessian
[Dxx,Dxy,Dyy] = my_Hessian2D(Image);

% note: abs(Lambda2)<abs(Lambda1)
[Lambda2,Lambda1,e2_x,e2_y]=eig2image(Dxx,Dxy,Dyy);
dimensions=size(Image);
lcs=zeros(dimensions);
%lcs=nan(dimensions);

%%{
%cond1=(Dx.*e2_x+Dy.*e2_y)<epsilon;

cond1=(Dx.*-e2_y+Dy.*e2_x)<epsilon;
cond2=(Lambda2<0);
%lcs(cond1)=1;
lcs_1=lcs;
lcs_1(cond1)=1;
%lcs_1(cond1&cond2)=1;
figure
pcolor(X, Y, lcs_1); colorbar; shading flat;
title('ones')
%}
%%{
%%
Hgx=Dxx.*Dx + Dxy.*Dy;
Hgy=Dxy.*Dx + Dyy.*Dy;

normHg=sqrt(Hgx.^2+Hgy.^2);

% note abs(Lambda2)<abs(Lambda1)

normLambda1g=abs(Lambda1).*sqrt((Dx).^2+(Dy).^2);
%normLambda1g=sqrt((Lambda1.*Dx).^2+(Lambda1.*Dy).^2);

cond1=(abs(normHg-normLambda1g)<epsilon);
cond2=(Lambda2<0);
lcs(cond1)=1;
%lcs(cond1&cond2)=1;
%}
fprintf(' zeros: %i\n ones:  %i\n nans:  %i\n',...
	length(find(lcs==0)),length(find(lcs==1)),length(find(isnan(lcs))))
%[X_grid,Y_grid]=meshgrid(X,Y);
figure
title('ones')
pcolor(X, Y, lcs); colorbar; shading flat;
end

function D = gradient2(F,option)
[k,l] = size(F);
D  = zeros(size(F),class(F)); 

switch lower(option)
case 'x'
    % Take forward differences on left and right edges
    D(1,:) = (F(2,:) - F(1,:));
    D(k,:) = (F(k,:) - F(k-1,:));
    % Take centered differences on interior points
    D(2:k-1,:) = (F(3:k,:)-F(1:k-2,:))/2;
case 'y'
    D(:,1) = (F(:,2) - F(:,1));
    D(:,l) = (F(:,l) - F(:,l-1));
    D(:,2:l-1) = (F(:,3:l)-F(:,1:l-2))/2;
otherwise
    disp('Unknown option')
end
        
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