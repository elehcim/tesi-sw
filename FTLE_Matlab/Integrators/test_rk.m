clear
mu=0.012277471;
ci=[1-mu-0.05 0 0.005 0.5290];
tspan=[0 2];
Y=RK4(@f,tspan,ci,mu);
figure
plot(Y(:,1),Y(:,2))
Y(end,1)