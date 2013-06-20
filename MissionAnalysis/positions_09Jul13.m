W_jup=274.05*pi/180; % rad
TA_jup=75.42*pi/180;
RG_jup= 7.667325090818796E+08; % km


W_earth=262.30*pi/180;
TA_earth=184.25*pi/180;
RG_earth=1.524498908899928E+08; % km


long_jup=W_jup+TA_jup;

long_earth=W_earth+TA_earth;

long_jup_d=mod(long_jup*180/pi,360)
long_earth_d=mod(long_earth*180/pi,360)

delta_angle=360-long_jup_d+long_earth_d

%% Synodic reference system
f_jup=TA_jup
L=RG_jup;

x_earth=(RG_earth*cosd(delta_angle))/L
y_earth=(RG_earth*sind(delta_angle))/L
