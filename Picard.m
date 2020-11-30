%% ODE45
p.explicit = true;
DeriveEOM;

steps = 10000;
t0 = 0;
tf = 10;

tspan = linspace(t0,tf,steps);
q0 = [2*ones(1,n) zeros(1,n)];

[t,q] = ode45(@(t,q) p.pendODE_e(t,q),tspan,q0);

%% Picard Operator
clear z;
z(:,:,1) = q0'*ones(1,steps);
i=1;
while 1
    i=i+1
    for t_ = 1:steps
        int_t0_t =  trapz(tspan(1:t_),p.pendODE_e(tspan(1:t_),z(:,1:t_,i-1)),2);
        z(:,t_,i) =  z(:,t_,1) + int_t0_t;
    end
    if norm(z(:,:,i) - z(:,:,i-1))<0.001
        break;
    end
end
max_iter = i;
%% Plot iterations
hold on
for j = 1:max_iter
plot(t,z(:,:,j)','Color',[(max_iter-j)/max_iter j/max_iter .1],'LineWidth',2)
end
% title('Picard Iterations for n=1','interpreter','latex')
ylabel('$[\theta, \dot{\theta}]$','interpreter','latex')
xlabel('$T$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',4)
%% Plot difference
hold on;
p1 = plot(t,z(:,:,end)','g','LineWidth',3);
p2 = plot(t,q,'k--','LineWidth',3);

% title('$Picard vs ode45 for n=1$','interpreter','latex')
ylabel('$[\theta, \dot{\theta}]$','interpreter','latex')
xlabel('$T$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',4)
legend([p1(1), p2(1)],{'$Picard$','ODE45'},'interpreter','latex','location','southwest')


%% Iterative Picard
steps = 101;
ti = 0.1;

z=q0';
% z = [];
for j=1:ceil(tf/ti)
clear z_i;
z_i(:,:,1) = z(:,end)*ones(1,steps);
disp(j);
i=1;
while 1
    i=i+1;
    for t_ = 1:steps
        int_t0_t =  trapz(tspan(1:t_),p.pendODE_e(tspan(1:t_),z_i(:,1:t_,i-1)),2);
        z_i(:,t_,i) =  z_i(:,t_,1) + int_t0_t;
    end
    if norm(z_i(:,:,i) - z_i(:,:,i-1))<0.001
        break;
    end
end
z = [z z_i(:,2:end,end)];
end
%% Plot difference
hold on;
p1 = plot(t,z(:,1:end-1)','g','LineWidth',3);
p2 = plot(t,q,'k--','LineWidth',3);

% title('$n=1$','interpreter','latex')
ylabel('$[\theta, \dot{\theta}]$','interpreter','latex')
xlabel('$T$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',4)
legend([p1(1), p2(1)],{'$Picard$','ODE45'},'interpreter','latex','location','northwest')


%% Appendix
%         int_t0_t = zeros(6,1);
%         for tau = 1:t
%             int_t0_t = int_t0_t + pendODE(tspan(tau),z(:,tau,i-1))*tf/steps;
%         end
% int_t0_t = int_t0_t + integral(@(tau) pendODE(tspan(tau),z(:,tau,i-1)),1,t,'ArrayValued',true,'Waypoints',1:t);
% int_t0_t = int_t0_t + integral(@(tau) pendODE(tau,interp1(tspan,z(:,:,i-1)',tau)'),0,tspan(t),'ArrayValued',true);
