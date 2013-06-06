clear all
tt=tic;

%% Computation parameters
n=200;
% Adesso la variabile indipendente è l'anomalia vera "ni"
% 1 periodo = 2*pi -> F = 2 è circa il 31% del periodo (1/pi di periodo)
F=2;
n_frame=10;
ni_span=linspace(0,2*pi,n_frame);

frame_time=zeros(n_frame);
folder_name=[folder, 'ER3BP/Figure/'];

%% Assigned initial conditions
ecc=0.04;
mu=0.1;
[xl1,yl1,xl2,yl2,xl3,yl3,xl4,yl4,xl5,yl5]=Lagr(mu);
C_L1=2*Omega(xl1,yl1,mu); %FIXME ci vuole /k?
E_0=-C_L1/2+0.03715;
Y_0=0;

%% Main cycle
for index=1:n_frame
	ni=ni_span(index);
	FTLE_Gawlick_ell
end
total_time=toc(tt);
disp('----Tempo impiegato----\n');
fprintf('Frame \t Tempo (min)\n');
for pippo=1:n_frame
	fprintf('%2i |\t %.2f\n',pippo,frame_time(pippo)/60);
end
fprintf('Tot\t |\t %.2f\n',total_time/60);