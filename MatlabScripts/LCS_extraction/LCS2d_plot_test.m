function [X_res,Y_res,lcs]=LCS2d_plot_test(Image, epsilon, X, Y, method, resample)
% lcs=LCS2d_plot_test(Image, epsilon,X,Y, method, resample)

%% resampleImage
if resample
	[X_res,Y_res,v1]=resample_image( Image, X, Y, 1000 );
	% 	figure
	% 	pcolor(v1);shading flat;colorbar; title('resampled image')
	Image=v1;
else
	X_res=X;
	Y_res=Y;
end

%% Compute Gradient
grad_comp='';
switch grad_comp
	case 'mine'
		Dx=gradient2(Image,'x');
		Dy=gradient2(Image,'y');
	otherwise
		[Dx,Dy]=gradient(Image);
end

%% Hessian
[Dxx,Dxy,Dyy] = my_Hessian2D(Image);

%% Implement method
dimensions=size(Image);
lcs=zeros(dimensions);
switch method
	case 'cross product'
		cross_product = (Dx.*Dy.*(Dxx-Dyy) + Dxy.*(Dy.^2-Dx.^2));
		cond=cross_product<epsilon;
		lcs(cond)=1;
		
	case 'curvature vector'
		[Lambda2,~,e2_x,e2_y]=eig2image(Dxx,Dxy,Dyy);
		cond1=(Dx.*e2_x+Dy.*e2_y)<epsilon;
		cond2=(Lambda2<0);
		lcs(cond1&cond2)=1;
		
	case 'eigen vector'
		% note abs(Lambda2)<abs(Lambda1)
		[Lambda2,Lambda1,~,~]=eig2image(Dxx,Dxy,Dyy);
		
		Hgx=Dxx.*Dx + Dxy.*Dy;
		Hgy=Dxy.*Dx + Dyy.*Dy;
		
		normHg=sqrt(Hgx.^2+Hgy.^2);
		normLambda1g=abs(Lambda1).*sqrt((Dx).^2+(Dy).^2);
		%normLambda1g=sqrt((Lambda1.*Dx).^2+(Lambda1.*Dy).^2);
		
		% FIXME - WHY turns out to be the opposite?????
		cond1=(abs(normHg-normLambda1g)<epsilon);
		cond2=(Lambda2<0);
		lcs(cond1&cond2)=1;
	otherwise
		error('please specify a method')
end
fprintf(' zeros: %i\n ones:  %i\n nans:  %i\n',...
	length(find(lcs==0)),length(find(lcs==1)),length(find(isnan(lcs))))
pcolor(X_res,Y_res,lcs); shading interp; colorbar;
title(sprintf('%s\n \\epsilon = %.0e',method,epsilon))

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

%{
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
%}