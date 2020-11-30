function p = params()
p = struct();
p.explicit= true;
p.redoEOM = true;
p.redoEnergy = true;

p.controls = true; % :)
p.controller = @controller_main;

%% Pendulum Parameters
p.n = 2;
p.g = 9.81;
p.m = 1*ones(1,p.n);
p.b = 0*ones(1,p.n);
p.l = 1*ones(1,p.n);
% p.l=rand(1,p.n);

%% Controller Parameters
p.kp = 100*ones(1,p.n);
p.kd = 30*ones(1,p.n);
p.q_des = pi*ones(1,p.n);
p.qd_des = zeros(1,p.n);

%% Simulation Parameters
p.tspan = linspace(0,0.14655,300);
% p.q0 = [pi/2*ones(p.n,1); 0*ones(p.n,1)];
% p.q0 = [pi*ones(p.n,1); 0.1*ones(p.n,1)];
% p.q0 = [3*pi/2; 3*pi/2; 0; 1];
p.q0 = [-1.3*pi/4; 0; 0; 6];

%% Animation Parameters
p.drawTail = true;

%% Lipschitz Connstants
p.epsilon = 1;
L0 = [p.g 40 76.8 119.8];
try
    p.L0 = L0(p.n);
catch E
    warning("L0 not yet generated for choice of n");
end

p.pendODE_e = eval(['@(t,q,p) pendODE_e_' num2str(p.n) '(t,q,p)']);
p.pendODE_i = eval(['@(t,q,dq,p) pendODE_i_' num2str(p.n) '(t,q,dq,p)']);
p.pendEnergy = eval(['@(q,p) energy_' num2str(p.n) '(q,p)-energy_' num2str(p.n) '(zeros(1,2*p.n),p)']);

p.filename_i = [];
p.filename_e = [];
p.energyName = [];
end