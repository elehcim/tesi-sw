% generate a set of tracers starting for x=-0.1922 and y=0 inside the LCS
% of energy e=-1.35 in the Sun-Jupiter system
clearvars
n=5;
tr.n_tracers=n;
tr.x=-0.1922*ones(1,n);
tr.y=zeros(1,n);
tr.vx=linspace(1.534,1.789,n); % these value are proper of the LCS
tr.e=-1.35*ones(1,n);
tr.t0=2.2;
tr.T=5;
tr.mu=9.537e-4;
tr.ecc=0.04839;
tr=complete_tracers(tr); % compute vy
disp(tr)