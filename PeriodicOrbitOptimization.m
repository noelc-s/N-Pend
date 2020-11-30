%%
p=params();

ord = 10;

p.q0 = [pi/4*ones(p.n-1,1);0; 0*ones(p.n,1)];
p.qf = p.q0;

u=getBezier(ord);
p.u = u;
p.controller = @controller;

tau = linspace(0,1);
p.tspan = linspace(0,1,100);

p.control_action = logical([1 0 0]);
p.actuation = sum(sum(p.control_action));

%%
options = optimoptions('fmincon','Display','iter','StepTolerance',1e-10,'FunctionTolerance',1e-10,'RelLineSrchBnd',1e-10,'UseParallel',true,'MaxFunctionEvaluations',1000);
tic
[coefs,val] = fmincon(@(coefs) getMetric(u,tau,coefs),[[3;zeros(p.n-(1+p.n-p.actuation),1)] p.control_action*1*ones(p.n,ord+1)],...
    [],[],[],[],[[0;zeros(p.n-(1+p.n-p.actuation),1)] -100*p.control_action*ones(p.n,ord+1)],[[100;zeros(p.n-(1+p.n-p.actuation),1)] 100*p.control_action*ones(p.n,ord+1)],...
    @(coefs) constraints(p,tau,coefs),options);
toc
%%
p.coefs = coefs;
p.tspan = linspace(0,coefs(1,1),200);
% p.tspan = linspace(0,coefs(1,1)*10,1000);
p.tau_end = coefs(1,1);
p.q0 = [pi/4*ones(p.n-1,1);0; 0*ones(p.n,1)];
% p.q0 = [0.5960    0.7559   -0.1068   -1.2674   -0.8829   -1.5935]';

p.controller = @(t,q,p) controller(t,q,p);
% p.controller = @(t,q,p) controller_main(mod(t,p.tau_end),q,p,'traj',q_des,'tspan',p.tspan_one');
% p.controller = @(t,q,p) controller(t,q,p)+ controller_main(mod(t,p.tau_end),q,p,'traj',q_des,'tspan',p.tspan_one');

[t,q] = solveODE(p);

plot(q)

%% Construct Poincare Section
p.q0 = [pi/4*ones(p.n-1,1);0; 0*ones(p.n,1)];
% p.q0 = [0.5960    0.7559   -0.1068   -1.2674   -0.8829   -1.5935]';
p.f_0 = p.pendODE_e(0,p.q0,p);
p.tspan = linspace(0,coefs(1,1)*2,200);

% [q0 - q]*p.f_0;

[x,y] = meshgrid((-0.5:0.1:0.5));
z = (p.f_0(4)*(p.q0(4)-x)+p.f_0(5)*(p.q0(5)-y))/p.f_0(6)+p.q0(6);

surf(x,y,z)
shading interp
hold on
quiver3(p.q0(4),p.q0(5),p.q0(6),p.f_0(4)/norm(p.f_0),p.f_0(5)/norm(p.f_0),p.f_0(6)/norm(p.f_0),'MaxHeadSize',2,'LineWidth',2)
plot3(q(:,4),q(:,5),q(:,6),'LineWidth',3)
scatter3(q(1,4),q(1,5),q(1,6))
axis equal
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'FontSize',17)
    set(gca,'linewidth',2)
    xlabel('$$\dot{\theta}_1$$','interpreter','latex')
    ylabel('$$\dot{\theta}_2$$','interpreter','latex')
    zlabel('$$\dot{\theta}_3$$','interpreter','latex')
%%
p.q0 = [pi/4*ones(p.n-1,1);0; 0*ones(p.n,1)];
% p.q0 = [0.5960    0.7559   -0.1068   -1.2674   -0.8829   -1.5935]';
% p.q0(2) = p.q0(2) + 0.1;
options = odeset('Events',@(t,q) PoincareStrike(t,q,p));

t2 = [0];
t = [0];
DI = [];
q_nom = p.q0';
q_sigma = p.q0';
for i = 1:2
    p.q0 = q_nom(end,:)';
    [t1,q1] = ode45(@(t,q) p.pendODE_e(t,q,p),p.tspan, p.q0,options);
    [t2,q2] = ode45(@(t,q) p.pendODE_e(t,q,p),p.tspan+t1(end),q1(end,:),options);
    t=[t;t(end)+[t1;t2]];
    DI = [DI;zeros(length(t1)+length(t2)-1,1); 1];
    q_nom=[q_nom;[q1;q2]];
    q_sigma = [q_sigma; q_nom(end,:)];
end
%% Plot Stable Orbit that shows up
figure('units','normalized','outerposition',[0 0 0.43 1])
plot3(q_nom(1:120,4),q_nom(1:120,5),q_nom(1:120,6),'b','LineWidth',2)
hold on
plot3(q_nom(end-500:end,4),q_nom(end-500:end,5),q_nom(end-500:end,6),'r','LineWidth',2)
axis equal
axis([min(q_nom(:,4)) max(q_nom(:,4)) min(q_nom(:,5)) max(q_nom(:,5)) min(q_nom(:,6)) max(q_nom(:,6))])
% axis off
set(gca, 'Color', 'none')
% plot(q_sigma)

%% Plot distance to orbit
T_ind = find(t>p.tau_end,1)-1;
O = q_nom(1:T_ind,:);
for i = 1:length(t)-T_ind
    i_ = i+T_ind-1;
    e(i,:) = q_nom(i_+1,:) - O(mod(i_-1,T_ind)+1,:);
    d(i) = norm(e(i,:));
    d_3d(i,:) = [d(i)*sin((i_-1)/T_ind*2*pi) d(i)*cos((i_-1)/T_ind*2*pi)];
    d_ZOH(i) = max(d(1:i)); 
end

[X,Y,Z] = cylinder(d_ZOH);
% plot3(log(t(1:length(t)-T_ind)),d_3d(:,1), d_3d(:,2),'LineWidth',2)
patch(log(t(1:length(t)-T_ind)),d_3d(:,1), d_3d(:,2),log(t(1:length(t)-T_ind)),'FaceColor','none','EdgeColor','interp')
hold on
surf(log(Z)+5.6,X,Y,log(Z),'EdgeAlpha',0,'FaceAlpha',0.2);



%% Save Video of Transversing

figure('units','normalized','outerposition',[0 0 0.43 1])

[x,y] = meshgrid(-0.5:0.1:0.5);
z = (p.f_0(4)*(p.q0(4)-x)+p.f_0(5)*(p.q0(5)-y))/p.f_0(6)+p.q0(6);

surf(x,y,z)
shading interp
hold on
% quiver3(p.q0(4),p.q0(5),p.q0(6),p.f_0(4)/norm(p.f_0),p.f_0(5)/norm(p.f_0),p.f_0(6)/norm(p.f_0),'MaxHeadSize',2,'LineWidth',2)
Xcoord = [];
Ycoord = [];
Zcoord = [];
Ccoord = [];
c = [];
clear M;
h=patch(q_nom(1,4),q_nom(1,5),q_nom(1,6),c(:),'FaceColor','none','EdgeColor','interp')
s = scatter3(q_nom(1,4),q_nom(1,5),q_nom(1,6),30,'r','filled');
axis equal
axis([min(q_nom(:,4)) max(q_nom(:,4)) min(q_nom(:,5)) max(q_nom(:,5)) min(q_nom(:,6)) max(q_nom(:,6))])
% axis off
set(gca, 'Color', 'none')
colormap(jet)
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'FontSize',17)
    set(gca,'linewidth',2)
    xlabel('$$\dot{\theta}_1$$','interpreter','latex')
    ylabel('$$\dot{\theta}_2$$','interpreter','latex')
    zlabel('$$\dot{\theta}_3$$','interpreter','latex')
for i = 1:2:size(q_nom,1)
    Xcoord = [Xcoord q_nom(i:i+1,4)'];
    Ycoord = [Ycoord q_nom(i:i+1,5)'];
    Zcoord = [Zcoord q_nom(i:i+1,6)'];
    Ccoord = [Ccoord 2*(i:i+1)/size(q_nom,1)-1];
    x = [Xcoord NaN];
    y = [Ycoord NaN];
    z = [Zcoord NaN];
    c = [Ccoord NaN];
    set(h,'XData',x,'YData',y,'ZData',z,'CData',c);
    set(s,'XData',Xcoord(end),'YData',Ycoord(end),'ZData',Zcoord(end));
    view((i/800)*50,30)
    M(ceil(i./2)) = getframe(gcf);

drawnow
end
%% P Sec closeup
figure('units','normalized','outerposition',[0 0 0.3 0.6])
[x,y] = meshgrid(-2:0.1:2);
z = (p.f_0(4)*(p.q0(4)-x)+p.f_0(5)*(p.q0(5)-y))/p.f_0(6)+p.q0(6);

surf(x,y,z)
shading interp
hold on
s = scatter3(q_sigma(1,4),q_sigma(1,5),q_sigma(1,6),30,'r','filled');
axis off
axis equal

Xcoord = [];
Ycoord = [];
Zcoord = [];
clear M;
j = 0;
for i = 1:size(q_nom,1)-1
    if (DI(i) == 1)
        j = j+1;
        Xcoord = [Xcoord q_sigma(j,4)];
        Ycoord = [Ycoord q_sigma(j,5)];
        Zcoord = [Zcoord q_sigma(j,6)];
        set(s,'XData',Xcoord,'YData',Ycoord,'ZData',Zcoord);
    end
    drawnow
        M(i) = getframe(gcf);

end


% Results: the main orbit is not stable, but after some transience, it
% levels off to a period 2 orbit. So cool

%% Energy?
for i =1:size(q_nom,1)
   E(i) = energy_3(q_nom(i,:),p); 
   Ed(i) = p.b*(q_nom(i,4:end).^2)' + q_nom(i,4:end)*p.controller(t(i),q_nom(i,:),p);
end

plot(t,-E,'linewidth',2);
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'FontSize',17)
    set(gca,'linewidth',2)
    xlabel('Time','interpreter','latex')
    ylabel('Energy','interpreter','latex')
    set(gca, 'Color', 'none')

%% Do perturbations
p.tspan = linspace(0,coefs(1,1),100);
[t,q] = solveODE(p);

q_nom = q;
dx = 0.01;

options = odeset('Events',@(t,q) PoincareStrike(t,q,p));

for i = 1:2*p.n
    
% p.q0 = [0.5960    0.7559   -0.1068   -1.2674   -0.8829   -1.5935]';
p.q0 = [pi/4*ones(p.n-1,1);0; 0*ones(p.n,1)];
    
t2 = [0];
t = [0];
DI = [];
p.q0(i) = p.q0(i)+dx;
q_per{i} = p.q0';
q_sigma{i} = p.q0';

for iter = 1:2
    p.q0 = q_per{i}(end,:)';
    [t1,q1] = ode45(@(t,q) p.pendODE_e(t,q,p),p.tspan, p.q0,options);
    [t2,q2] = ode45(@(t,q) p.pendODE_e(t,q,p),p.tspan+t1(end),q1(end,:),options);
    t=[t;t(end)+[t1;t2]];
    DI = [DI;zeros(length(t1)+length(t2)-1,1); 1];
    q_per{i}=[q_per{i};[q1;q2]];
    q_sigma{i} = [q_sigma{i}; q_per{i}(end,:)];
end
end

%%
v = 2;
p.q0 = [0.5960    0.7559   -0.1068   -1.2674   -0.8829   -1.5935]';
% p.q0 = [pi/4*ones(p.n-1,1);0; 0*ones(p.n,1)];

for i = 1:2*p.n
    for j = 1:2*p.n
        DP(i,j) =  (q_sigma{i}(v,j) - q_sigma{i}(1,j));
    end
end

plot(q_sigma{1},'r')
hold on;
plot(q_sigma{2},'g')
plot(q_sigma{3},'b')
plot(q_sigma{4},'k')
plot(q_sigma{5},'y')
plot(q_sigma{6},'c')

%%
function [value,isterminal,direction] = PoincareStrike(t,q,p)
x0 = p.q0';
x = q';
value = (x0-x)*p.f_0;
isterminal = 1;   % Stop the integration
direction = -1;   % Negative direction only
end

function metric = getMetric(u,tau,coefs)
metric = 0;
for i = 1:size(coefs,1)
    metric = metric + norm(u(tau,coefs(i,2:end)));
end
end

function [c,ceq] = constraints(p,tau,coefs)
p.coefs = coefs;
p.tspan = linspace(0,coefs(1,1),100);
[t,q] = solveODE(p);
c = -1;
ceq = q(end,:) - p.qf';
end

function u=controller(t,q,p)
p.tau_end = p.coefs(1,1);
u = [];
tau = mod(t/p.tau_end,1);
for i = 1:p.actuation
    u = [u p.u(tau,p.coefs(i,2:end))];
end
for i = 1:(p.n-p.actuation)
    u = [u zeros(size(u,1),1)];
end
u=u';
% u(:,2:end) = 0;
end
%%
function B = getBezier(n)
B = @(t,coefs) (1-t).^n*coefs(1);
for i = 1:n
    tmp = char(B);
    eval(['B = @(t,coefs)' tmp(11:end) '+ nchoosek(' num2str(n) ',' num2str(i) ').*(1-t).^(' num2str(n) '-' num2str(i) ').*t.^' num2str(i) '*coefs(' num2str(i) '+1);']);
end
end