Driver
%% IC

Delta = p.q0;
Delta = Delta*2;

q0 = Delta;

n = p.n;
tspan = p.tspan;

dq0 = zeros(2*n,1);

[y0_new,yp0_new] = decic(@(t,q,dq) p.pendODE_i(t,q,dq),0,q0,[],dq0,[]);

options = odeset('RelTol',.1,'AbsTol',.1);

[t_,q_] = ode15i(p.pendODE_i,tspan,y0_new,yp0_new,options);

subplot(2,1,1);
hold on;
% plot(t_,q_(:,1),'LineWidth',2,'Color',[1 0 0]);
% plot(t_,q_(:,2),'LineWidth',2,'Color',[.5 0 0]);
% plot(t,q(:,1),'--','LineWidth',2,'Color',[0 0 1]);
% plot(t,q(:,2),'--','LineWidth',2,'Color',[0 0 .5]);
p1 = plot(t,q,'r','LineWidth',2);
p2 = plot(t_,q_,'b','LineWidth',2);


ylabel('$[\theta, \dot{\theta}]$','interpreter','latex')
xlabel('$T$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',4)
legend([p1(1),p2(1)],{'$q$','$\tilde{q}$'},'interpreter','latex','location','northeast')
subplot(2,1,2);

hold on;
plot(t_,sqrt(sum((q_-q).^2,2)),'LineWidth',2);
plot(t_,norm(p.q0-q0)*exp(p.L0*t),'LineWidth',2);
ylabel('$\|\cdot\|$','interpreter','latex')
xlabel('$T$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',4)
legend({'$\|x(t)-z(t)\|$','$\|x(0)-z(0)\|e^{Lt}$'},'interpreter','latex','location','northeast')

%% Parameters
% p.n = 1;
p.g = 19.81;
p.m = .1;
p.b = 0;
p.l = .1;
p.filename_i = 'tmp_i_';
p.filename_e = 'tmp_e_';
p.explicit = true;

DeriveEOM(p);

p.pendODE_i = eval(['@(t,q,dq) ' p.filename_i num2str(p.n) '(t,q,dq)']);
p.pendODE_e = eval(['@(t,q) ' p.filename_e num2str(p.n) '(t,q)']);

[t_,q_] = solveODE(p);
g = @(q) -norm(p.pendODE_e(t,q) - pendODE_e_1(t,q));

[x,mu] = fmincon(g,p.q0,[],[],[],[],-1*ones(1,2*p.n),1*ones(1,2*p.n));
% [x,mu] = fmincon(g,p.q0);

subplot(2,1,1);
hold on;
p1 = plot(t_,q_,'r','LineWidth',2);
p2 = plot(t,q,'b','LineWidth',2);
ylabel('$[\theta, \dot{\theta}]$','interpreter','latex')
xlabel('$T$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',4)
legend([p1(1),p2(1)],{'$q$','$\tilde{q}$'},'interpreter','latex','location','northeast')

subplot(2,1,2);
hold on;
plot(t_,sqrt(sum((q_-q).^2,2)),'LineWidth',2);
plot(t_,-mu*t.*exp(p.L0*t),'LineWidth',2);
ylabel('$\|\cdot\|$','interpreter','latex')
xlabel('$T$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',4)
legend({'$\|x(t)-z(t)\|$','$\mu te^{Lt}$'},'interpreter','latex','location','northeast')