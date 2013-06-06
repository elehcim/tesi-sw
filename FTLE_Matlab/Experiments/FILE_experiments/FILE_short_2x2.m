clear all
tic
%% Computation parameters
n=100;
T=30;
N=5;
%% Initial conditions ranges

% Assigned initial conditions
mu=0.012151; %Earth-Moon
[xl1,yl1,xl2,yl2,xl3,yl3,xl4,yl4,xl5,yl5]=Lagr(mu);
C_L2=2*Potential(xl2,yl2,mu);

x_0_min=0.1;
x_0_max=0.6;

y_0=0;

vx_0_min=-0.55;
vx_0_max=-0.25;

x_0=linspace(x_0_min, x_0_max,n);
vx_0=linspace(vx_0_min, vx_0_max,n);

dx=x_0(2)-x_0(1);
dvx=vx_0(2)-vx_0(1);
% Compute initial conditions
vy_0=zeros(n,n);
for i=1:n
    for j=1:n
        vy_0(i,j)=sqrt(2*Potential(x_0(i),y_0,mu)-C_L2-vx_0(j)^2);
    end
end

%% Look for phisically meaningful points
filter=zeros(n,n);  %0=meaningless point 1=meaningful point
E_k=-ones(n,n);
for i=1:n
    for j=1:n
        if isreal(vy_0(i,j))
        E_k(i,j)=0.5*(vx_0(j)^2+vy_0(i,j)^2);
            if E_k(i,j)>=0
                filter(i,j)=1;
            end
        end
    end
end

%% Integration
x_T=NaN*ones(n,n);
y_T=NaN*ones(n,n);
vx_T=NaN*ones(n,n);
vy_T=NaN*ones(n,n);

%options=odeset('AbsTol',1e-12,'RelTol',1e-12,'Events',@cross_y); % From JD James Mireles
options=odeset('Events',@cross_y); % From JD James Mireles
h=waitbar(0,'','Name','Integration in progress, please wait!');

for i=1:n
    waitbar(i/n,h,sprintf('Computing i=%i',i));
    c1=x_0(i);
    parfor j=1:n
        if filter(i,j)
            [t,Y,te,ye,ie]=ode45(@f,[0 T],[c1; y_0; vx_0(j); vy_0(i,j)],options,mu);
            x_T(i,j)=ye(N+1,1);
            vx_T(i,j)=ye(N+1,3);
            y_T(i,j)=ye(N+1,2);
            vy_T(i,j)=ye(N+1,4);
        end
    end
end
close(h);
%% FTLE Computation
filter_ftle=filter;

% Border exclusion
filter_ftle(1,:)=0;
filter_ftle(n,:)=0;
filter_ftle(:,1)=0;
filter_ftle(:,n)=0;

% Exclusion of meningless point and their neighbours
for i=2:(n-1)
    for j=2:(n-1)
        if filter(i,j)==0
            filter_ftle(i+1,j)=0;
            filter_ftle(i-1,j)=0;
            filter_ftle(i,j+1)=0;
            filter_ftle(i,j-1)=0;
        end
    end
end

% FTLE Computation

file=NaN*ones(n,n);

for i=1:n
    for j=1:n
        if filter_ftle(i,j)
            dphi(1,1)=(x_T(i+1,j)-x_T(i-1,j))/dx;
			
			dphi(1,2)=(x_T(i,j+1)-x_T(i,j-1))/dvx;
	
			dphi(2,1)=(vx_T(i+1,j)-vx_T(i-1,j))/dx;
            
			dphi(2,2)=(vx_T(i,j+1)-vx_T(i,j-1))/dvx;
            
            file(i,j)=(1/abs(N))*log(norm(dphi));
        end
    end
end

%% Results' plot
figure
pcolor(x_0,vx_0,file')
shading flat
t_tot=toc
nome=['var_xvx_', 'ode00', '_n',num2str(n),'_FILE_short'];
save(nome)