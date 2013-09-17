function f = calc_t0( nu,GM,ecc,omega, a_jup )
delta_t=nu/omega;
eps=0.001;
n=sqrt(GM/(a_jup^3));
M=n*(delta_t);
E=M;
E_new=E;
i=0;
while i<10
	E=E_new;
	E_new=E+(M-E+ecc*sin(E))/(1-ecc*cos(E));
	i=i+1;
	%fprintf('E    =   %.10e \n',E);
	%fprintf('E_new=   %.10e \n',E_new);
	if E_new-E<eps
		break
	end
end
	f=2*atan2(sqrt(1+ecc)*sin(E_new/2),sqrt(1-ecc)*cos(E_new/2));