clear;close all;clc;
addpath(genpath(pwd));
p = params();

% DeriveEOM(p);
[t,q] = solveODE(p);
animatePend
%%

t_1=q(:,1);
td_1=q(:,3);
t_2=q(:,2);
td_2=q(:,4);
fk = [sin(t_1)+sin(t_2) -cos(t_1)-cos(t_2)];
h = 1-sin(t_1)-sin(t_2);
Jh = [-cos(t_1) -cos(t_2)];

alpha = 10;
alpha_e = 10;

for i = 1:length(t)
t_1=q(i,1);
td_1=q(i,3);
t_2=q(i,2);
td_2=q(i,4);
D = [(p.m(1)+p.m(2))*p.l(1)^2+p.m(2)*p.l(2)^2+2*p.m(2)*p.l(1)*p.l(2)*cos(t_2-t_1) p.m(2)*p.l(2)^2+p.m(2)*p.l(1)*p.l(2)*cos(t_2-t_1);p.m(2)*p.l(2)^2+p.m(2)*p.l(1)*p.l(2)*cos(t_2-t_1) p.m(2)*p.l(2)^2];
C = [0 -p.m(2)*p.l(1)*p.l(2)*(2*td_1+(td_2-td_1))*sin(td_2-td_1); p.m(2)*p.l(1)*p.l(2)*td_1*sin(t_2-t_1) 0];
G = p.g*[(p.m(1)+p.m(2))*p.l(1)*sin(t_1)+p.m(2)*p.l(2)*sin(t_2); p.m(2)*p.l(2)*sin(t_2)];
    % Is this right?
    % https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-832-underactuated-robotics-spring-2009/readings/MIT6_832s09_read_appA.pdf
h_D(i) = -0.5*[td_1 td_2]*D*[td_1; td_2]+alpha_e*h(i);
B = [0;1];

u_des(i)=-100*(t_2-pi)-10*(td_2);

Q = 1;
R = -u_des(i);
A(i) = [td_1 td_2]*B;
b = alpha*h_D(i) + alpha_e*Jh(i,:)*[td_1; td_2] + G'*[td_1; td_2];

u(i)=quadprog(Q,R,A(i),b);

val(i) = A(i)*u(i);
ub(i) = b;

end

%%
subplot(211); hold on
plot(t,h_D,'Linewidth',3)
plot(t,A,'Linewidth',3)
yline(0)
axis([0 t(end) -10 15])
xlabel('Time','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)
legend({'$h_{D}$','$\|\dot{q}^TB\|$'},'interpreter','latex');
subplot(212); hold on;
plot(t,-val,'Linewidth',10)
hold on
plot(t,-A.*u_des,'Linewidth',3)
plot(t,-ub,'Linewidth',3)
xlabel('Time','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)
legend({'$-\dot{q}^TBu$','$\dot{q}^TBu_{des}$','$LB$'},'interpreter','latex');