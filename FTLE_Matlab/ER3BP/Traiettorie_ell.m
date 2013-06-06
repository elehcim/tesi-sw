clear all
%mu=0.012151; %Earth-Moon
mu=0.012277471 %Earth-Moon after Mireles
[xl1,yl1,xl2,yl2,xl3,yl3,xl4,yl4,xl5,yl5] = Lagr(mu);
C=3.17;
%C=2*Potential(xl1,yl1,mu);
x_0=1;
y_0=0;
vx_0=0;
vy_0=sqrt(-C-vx_0^2+2*Potential(x_0,y_0,mu));
%vy_0=0;
F=2; % significa 2 periodi
T=2; % significa 2 periodi
C_star=2*Potential(x_0,y_0,mu)-(vx_0^2+vy_0^2);
%C=-(vx_0^2+vy_0^2)+2*Omega(x_0,y_0,mu);
E=-C/2;
ecc=0;

options=odeset('AbsTol',1e-22,'RelTol',1e-13); % From JD James Mireles


%Integro l'orbita con modello ellittico
[t_ell,Y_ell]=ode45(@f_ell,[0 F],[x_0; y_0; vx_0; vy_0],options,mu,ecc);
x_ell=Y_ell(:,1);
y_ell=Y_ell(:,2);
vx_ell=Y_ell(:,3);
vy_ell=Y_ell(:,4);

% Precisione orbita
delta_E_ell=abs(Energy(x_ell,y_ell,vx_ell,vy_ell,mu)-E);

% figure
% title('Precisione orbita')
% plot(delta_E0)

% Integro l'orbita con modello circolare
[t_circ,Y_circ]=ode45(@f,[0 T],[x_0; y_0; vx_0; vy_0],options,mu);
x_circ=Y_circ(:,1);
y_circ=Y_circ(:,2);
vx_circ=Y_circ(:,3);
vy_circ=Y_circ(:,4);

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
figure
hold all

% Regioni di Hill
contour(x,y,z,[C/2,C/2])
text(-2,-2,sprintf('C=%.2f',C))

% Attrattori
plot(-mu,0,'ok')
plot(1-mu,0,'ok')

% Punto iniziale e finale
plot(x_ell(1),y_ell(1),'sg')
plot(x_ell(end),y_ell(end),'sr')

h=plot(x_ell,y_ell,'b',x_circ,y_circ,'-.r');
legend(h,'ell','circ')