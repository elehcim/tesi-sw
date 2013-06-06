function dv=deiperbolize(vx,vy,mu,ni,ecc,a)
% Take absolute velocities and returns delta v required to make the orbit
% closed around the primary (Take as a reference Cassini Missions
vx_planet=sqrt(mu/(a*(1-ecc^2)))*(-sin(ni));
vy_planet=sqrt(mu/(a*(1-ecc^2)))*(ecc+cos(ni));
vx_rel=vx-vx_planet;
vy_rel=vy-vy_planet;
v_excess=(vx_rel^2+vy_rel^2)/2; %FIXME sqrt
dv=v_excess;
end
% Impulse at periapse similar to Cassini Mission
% rp=a*(1-e)