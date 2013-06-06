function [x,y,vx,vy] = RKF45_gpu(t0, T, x_0,y_0,vx_0,vy_0, mu )
%RK4 4th-order Runge Kutta integrator for GPU computation
%   This RK4 integrator is implemented to run on GPU. This limits the usage
%   of arrays in the function. Every input and output parameter is a
%   scalar. The function f_gpu is the vector field with scalar i/o
%   parameters.

eps_max=1e-7;
eps_min=1e-5;
eps_abs=1e-6;
h_max=100;
h_min=1e-7;
h=h_max;
t=t0;
x=x_0; y=y_0; vx=vx_0; vy=vy_0;
if ~isnan(vy)
	while t<T
		[k11, k12, k13, k14]=f_gpu(t,x,y,vx,vy,mu);
		c1=h/4;
		[k21, k22, k23, k24]=f_gpu(t+h/4, x+k11*c1, y+k12*c1, vx+k13*c1, vy+k14*c1,mu);
		c1=3/32*h; c2=9/32*h;
		[k31, k32, k33, k34]=f_gpu(t+h*3/8, x+k11*c1+k21*c2, y+k12*c1+k22*c2, vx+k13*c1+k23*c2, vy+k14*c1+k24*c2,mu);
		c1=1932/2197*h; c2=-7200/2197*h; c3=7296/2197*h;
		[k41, k42, k43, k44]=f_gpu(t+12/13*h,x+k11*c1+k21*c2+k31*c3, y+k12*c1+k22*c2+k32*c3, vx+k13*c1+k23*c2+k33*c3, vy+k14*c1+k24*c2+k34*c3,mu);
		c1=439/216*h; c2=-8*h; c3=3680/513*h; c4=-845/4104*h;
		[k51, k52, k53, k54]=f_gpu(t+h, x+k11*c1+k21*c2+k31*c3+k41*c4, y+k12*c1+k22*c2+k32*c3+k42*c4, vx+k13*c1+k23*c2+k33*c3+k43*c4, vy+k14*c1+k24*c2+k34*c3+k44*c4,mu);
		c1=-8/27*h; c2=2*h; c3=-3544/2565*h; c4=1859/4104*h; c5=-11/40*h;
		[k61, k62, k63, k64]=f_gpu(t+h/2, x+k11*c1+k21*c2+k31*c3+k41*c4+k51*c5, y+k12*c1+k22*c2+k32*c3+k42*c4+k52*c5, vx+k13*c1+k23*c2+k33*c3+k43*c4+k53*c5, vy+k14*c1+k24*c2+k34*c3+k44*c4+k54*c5,mu);
		
		x4=x+(25/216*k11+1408/2565*k31+2197/4104*k41-1/5*k51)*h;
		y4=y+(25/216*k12+1408/2565*k32+2197/4104*k42-1/5*k52)*h;
		vx4=vx+(25/216*k13+1408/2565*k33+2197/4104*k43-1/5*k53)*h;
		vy4=vy+(25/216*k14+1408/2565*k34+2197/4104*k44-1/5*k54)*h;
		
		
		x5=x+(16/135*k11+6656/12825*k31+28561/56430*k41-9/50*k51+2/55*k61)*h;
		y5=x+(16/135*k12+6656/12825*k32+28561/56430*k42-9/50*k52+2/55*k62)*h;
		vx5=x+(16/135*k13+6656/12825*k33+28561/56430*k43-9/50*k53+2/55*k63)*h;
		vy5=x+(16/135*k14+6656/12825*k34+28561/56430*k44-9/50*k54+2/55*k64)*h;
		
		eps=abs(x5-x4);
		if abs(y5-y4)>eps
			eps=abs(y5-y4);
		end
		if abs(vx5-vx4)>eps
			eps=abs(vx5-vx4);
		end
		if abs(vy5-vy4)>eps
			eps=abs(vy5-vy4);
		end
		
		if (eps<=eps_max) && (eps>=eps_min) % Error within imposed limits
			x=x4;
			y=y4;
			vx=vx4;
			vy=vy4;
			t=t+h;
		elseif eps>eps_max % Error above upper limit -> discard result and decrease step
			s=(eps_abs*h/(2*eps))^(1/4);
			if (s*h)>=h_min % New step above lower limit
				h=s*h;
			else % New step below lower limit -> cannot reduce step -> keep result and min step
				x=x4;
				y=y4;
				vx=vx4;
				vy=vy4;
				h=h_min;
				t=t+h;
			end
		else % Error below lower limit -> keep result and increase step
			x=x4;
			y=y4;
			vx=vx4;
			vy=vy4;
			s=(eps_abs*h/(2*eps))^(1/4);
			if (s*h)<=h_max % New step below upper limit -> adjust step
				h=s*h;
			else
				h=h_max; % New step above upper limit -> cannot increase step
			end
			t=t+h;
		end
	end
end
end
