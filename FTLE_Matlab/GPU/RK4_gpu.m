function [x,y,vx,vy] = RK4_gpu(t0, T, x_0,y_0,vx_0,vy_0, mu )
%RK4 4th-order Runge Kutta integrator for GPU computation
%   This RK4 integrator is implemented to run on GPU. This limits the usage
%   of arrays in the function. Every input and output parameter is a
%   scalar. The function f_gpu is the vector field with scalar i/o
%   parameters.
h=1e-3;
t=t0;
i=1;
x=x_0; y=y_0; vx=vx_0; vy=vy_0;
if ~isnan(vy)
while t<T
	[k11, k12, k13, k14]=f_gpu(t,x,y,vx,vy,mu);
	[k21, k22, k23, k24]=f_gpu(t+h/2,x+k11*h/2,y+k12*h/2,vx+k13*h/2,vy+k14*h/2,mu);
	[k31, k32, k33, k34]=f_gpu(t+h/2,x+k21*h/2,y+k22*h/2,vx+k23*h/2,vy+k24*h/2,mu);
	[k41, k42, k43, k44]=f_gpu(t+h,x+h*k31,y+h*k32,vx+h*k33,vy+h*k34,mu);
	x=x+(h/6*(k11+2*k21+2*k31+k41));
    y=y+(h/6*(k12+2*k22+2*k32+k42));
    vx=vx+(h/6*(k13+2*k23+2*k33+k43));
    vy=vy+(h/6*(k14+2*k24+2*k34+k44));
	t=t+h;
	i=i+1;
end 
end
end
