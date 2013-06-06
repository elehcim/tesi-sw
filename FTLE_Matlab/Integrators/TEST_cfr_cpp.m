clear
%mu=0.012151; %Earth-Moon
mu=0.012277471; %Earth-Moon after Mireles
x_0=1-mu-0.05;
y_0=0;
vx_0=0.005;
vy_0=0.5290;
T=2;

options=odeset('AbsTol',1e-22,'RelTol',1e-13); % From JD James Mireles

% Integro con RK4
tic
[Y1]=RK4(@f,[0 T],[x_0; y_0; vx_0; vy_0],mu);
t_RK4=toc
x1=Y1(:,1);
y1=Y1(:,2);
vx1=Y1(:,3);
vy1=Y1(:,4);
l1=length(Y1);


%Calcolo regioni di Hill
points=500;
bb=3; % Bounding box
x=linspace(-bb,bb,points);
y=linspace(-bb,bb,points);
[x,y]=meshgrid(x,y);
z=(Potential(x,y,mu));
% figure
% surfc(x,y,z,'Edgecolor','none')

%% Plot orbita
%figure
hold on
%contour(x,y,z,[C/2,C/2])
plot(x1,y1,'dr')
%text(-2,-2,sprintf('C=%.2f',C))
%plotto attrattori
plot(-mu,0,'ok')
plot(1-mu,0,'ok')
% Plot punto iniziale e finale
plot(x1(1),y1(1),'sg')