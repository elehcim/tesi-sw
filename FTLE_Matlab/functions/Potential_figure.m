mu=0.3;
%[x,y]=meshgrid(-1:0.01:1);
%z=Potential(x,y,mu);

points=500;
bb=3; % Bounding box
x=linspace(-bb,bb,points);
y=linspace(-bb,bb,points);
[x,y]=meshgrid(x,y);
z=(Potential(x,y,mu));
for i=1:500
	for j=1:500
		if z(i,j)>2
			z(i,j)=nan;
		end
	end
end
[xl1,yl1,xl2,yl2,xl3,yl3,xl4,yl4,xl5,yl5] = Lagr(mu);
P1=Potential(xl1,yl1,mu);
P2=Potential(xl2,yl2,mu);
P3=Potential(xl3,yl3,mu);
P4=Potential(xl4,yl4,mu);
P5=Potential(xl5,yl5,mu);
%mesh(x,y,-z)
hold on
surf(x,y,-z,'Edgecolor','none')
shading interp
plot3(xl1,yl1,-P1,'ro','MarkerSize',3,'MarkerFaceColor','r')
plot3(xl2,yl2,-P2,'ro','MarkerSize',3,'MarkerFaceColor','r')
plot3(xl3,yl3,-P3,'ro','MarkerSize',3,'MarkerFaceColor','r')
plot3(xl4,yl4,-P4,'ro','MarkerSize',3,'MarkerFaceColor','r')
plot3(xl5,yl5,-P5,'ro','MarkerSize',3,'MarkerFaceColor','r')