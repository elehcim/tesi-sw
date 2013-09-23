clear
dist_Sun_Jup = 778547200;	% km
r_earth_orbit = 149600000;	% km
R = r_earth_orbit/dist_Sun_Jup;
a_jup=778412027; %km
GM_jup=126711995; %km^3/s^2
GM_sun=132712439935; %km^3/s^2
GM=GM_jup+GM_sun;
ecc=0.04839;
v_E=29.783; % km/s
omega=2*pi/(365*24*3600); % earth omega in rad/s

%% Dati al perielio di Giove 2011-03-17T19:00
W_jup_p=274.02;
W_earth_p=288.29;
nu_earth_p=76.15;
lambda_earth_p=W_earth_p+nu_earth_p;

alpha_p=lambda_earth_p-W_jup_p

%% Dati al 9 Jul
W_jup_9Jul=274.05; % deg
nu_jup_9Jul=75.42;

W_earth_9Jul=262.30;
nu_earth_9_Jul=184.25;

lambda_jup_9Jul=W_jup_9Jul+nu_jup_9Jul;
lambda_earth_9Jup=W_earth_9Jul+nu_earth_9_Jul;

alpha_9Jul=lambda_earth_9Jup-lambda_jup_9Jul % Torna con il delta_angle dello script positions_09Jul13

nu=(alpha_9Jul+nu_jup_9Jul)*pi/180; % rad
giri=2; % infatti dal 2011 al 2013 sono 2 anni!!!
fprintf('nu della terra corrispondente al 9 luglio 2013 è %.2f rad\n',nu+giri*2*pi)
nu_jup_9Jul=75.42
calc_t0(nu+giri*2*pi,GM,ecc,omega, a_jup)*180/pi