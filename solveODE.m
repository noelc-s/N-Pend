function [t,q] = solveODE(p)
n = p.n;
tspan = p.tspan;
q0 = p.q0;

[t,q] = ode45(@(t,q) p.pendODE_e(t,q,p),tspan, q0);

% dq0 = zeros(2*n,1);
% y0_new = q0;
% yp0_new = dq0;

% options = odeset('RelTol',.001,'AbsTol',.001);
% [t,q] = ode15i(@(t,q,dq) p.pendODE_i(t,q,dq,p),tspan,y0_new,yp0_new,options);

% disp('System Simulated with ODE45');
end