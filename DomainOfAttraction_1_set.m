clear;close all;clc;
addpath(genpath(pwd));
p = params();
DeriveEOM(p);


[~,V,dV] = DeriveEnergy(p);
[t,q] = solveODE(p);

vars = {'q'};

for i = 1:2*p.n
    for var=vars
        eval(['syms ' var{:} '_' num2str(i)]);
    end
end

for i = 1:p.n
    V = subs(V,['b_' num2str(i)],p.b);
    V = subs(V,['m_' num2str(i)],p.m);
    V = subs(V,['l_' num2str(i)],p.l);
end

eq = solve(V==0,[q_1,q_2]);
%%
r = 2;
steps = 100;
theta = linspace(0,2*pi,steps);
ic = [r*cos(theta); r*sin(theta)]';
p.tspan = linspace(0,-20,1000);
for i = 1:steps
    p.q0 = ic(i,:);
    [t_tmp,q_tmp] = solveODE(p);
    Q(:,:,i) = q_tmp;
end
%%
hold on;
for i = 1:steps
plot(Q(:,1,i),Q(:,2,i))
drawnow;
% pause(0.05)
end
%%
hold on;
scatter(pi,0,200,'r','filled')
scatter(0*pi,0,200,'g','filled')
scatter(2*pi,0,200,'g','filled')
for i = 1:steps
% m1 = mod(Q(:,1,i)-pi,2*pi);
m1 = mod(Q(:,1,i),2*pi);
m2 = Q(:,2,i);
eps = .3;
m1(m1<eps)=NaN;
plot(m1,m2)
drawnow;
% pause(0.05)
end
axis([-0.2 2*pi+0.2 -18 18])

%%
hold on;
scatter3(0,0,1,200,'r','filled')
scatter3(0,0,-1,200,'g','filled')
for i = 1:steps
plot3(Q(:,2,i),sin(Q(:,1,i)),-cos(Q(:,1,i)))
drawnow;
view(50,30)
pause(0.05)
end