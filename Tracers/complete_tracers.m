function new_tr=complete_tracers(tr)
% Complete the tracers' missing field

%% see if some field is missing
vars={'x','y','vx','vy','e'};
a=isfield(tr,vars);
if sum(a)<4
	error('I can''t complete the tracer structure')
elseif sum(a)==4
	ind=find(a==0);
	tr.(vars{ind})=nan(1,tr.n_tracers);
end

%% Complete tracers' structure
first_tracer=[tr.x(1) tr.y(1) tr.vx(1) tr.vy(1) tr.e(1)];
if any(isnan(first_tracer))
	% see if the field is NaN
	index=find(isnan(first_tracer));
else
	warning('Nothing to complete')
	new_tr=tr;
	return
end
switch index
	case 3 % vx is missing
		tr.vx=-sign(tr.y).*sqrt(2*Omega(tr.x,tr.y,tr.mu)/(1+tr.ecc*cos(tr.t0))+2*tr.e-tr.vy.^2);
	case 4 % vy is missing
		tr.vy=sign(tr.x).*sqrt(2*Omega(tr.x,tr.y,tr.mu)/(1+tr.ecc*cos(tr.t0))+2*tr.e-tr.vx.^2);
	case 5 % e is missing
		tr.e=0.5*(tr.vx.^2+tr.vy.^2)-(Omega(tr.x,tr.y,tr.mu)/(1+tr.ecc*cos(tr.t0)));
end
new_tr=tr;