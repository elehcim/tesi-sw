mu=9.537e-4;
mu=0.1;
[xl1,yl1,xl2,yl2,xl3,yl3,xl4,yl4,xl5,yl5] = Lagr(mu);
figure
hold on
grid on
plot(0,0,'.k')
plot(-mu,0,'ok','MarkerFaceColor','k','Markersize',30)
plot(1-mu,0,'ok','MarkerFaceColor','k','Markersize',20)
plot(xl1,yl1,'rx','MarkerSize',30,'LineWidth',3)
plot(xl2,yl2,'rx','MarkerSize',30,'LineWidth',3)
plot(xl3,yl3,'rx','MarkerSize',30,'LineWidth',3)
plot(xl4,yl4,'r+','MarkerSize',30,'LineWidth',3)
plot(xl5,yl5,'r+','MarkerSize',30,'LineWidth',3)
plot([-1.5,1.5],[0,0],'k')
plot([0,0],[-1,1],'k')