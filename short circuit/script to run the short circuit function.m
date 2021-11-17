tspan=[0 300];
inicond=[0.15;0.033;0.75;0.04;1;293.15];

[x,T]=ode23s(@(x,T) libfunc(x,T) , tspan, inicond);
plot(x,T(:,6));


