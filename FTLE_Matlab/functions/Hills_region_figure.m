mu=0.3;
points=500;
bb=3; % Bounding box
x=linspace(-bb,bb,points);
y=linspace(-bb,bb,points);
[x,y]=meshgrid(x,y);
z=(Potential(x,y,mu));
[xl1,yl1,xl2,yl2,xl3,yl3,xl4,yl4,xl5,yl5] = Lagr(mu);
P1=Potential(xl1,yl1,mu);
P2=Potential(xl2,yl2,mu);
P3=Potential(xl3,yl3,mu);
P4=Potential(xl4,yl4,mu);
P5=Potential(xl5,yl5,mu);
figure
hold on
axis equal
axis ([-2 2 -2 2])
set(gca,'xtick',[])
set(gca,'ytick',[])
box on
plot(xl1,yl1,'ro','MarkerSize',5,'MarkerFaceColor','r','LineWidth',3)
plot(xl2,yl2,'ro','MarkerSize',5,'MarkerFaceColor','r','LineWidth',3)
plot(xl3,yl3,'ro','MarkerSize',5,'MarkerFaceColor','r','LineWidth',3)
plot(xl4,yl4,'ro','MarkerSize',5,'MarkerFaceColor','r','LineWidth',3)
plot(xl5,yl5,'ro','MarkerSize',5,'MarkerFaceColor','r','LineWidth',3)
plot(-mu,0,'ko','MarkerSize',20,'MarkerFaceColor','k')
plot(1-mu,0,'ko','MarkerSize',10,'MarkerFaceColor','k')
contour(x,y,z,[(P4)-0.1],'LineWidth',6,'Color',[0,0.188, 0.668])