tspan=[0 10000];
inicond=[0.15;0.033;0.75;0.04;300];

[x,T]=ode23s(@(x,T) cldfunc(x,T) , tspan, inicond);
plot(x,T(:,5));