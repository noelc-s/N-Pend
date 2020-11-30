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
% 
% for i = 1:p.n
%     V = subs(V,['b_' num2str(i)],p.b);
%     V = subs(V,['m_' num2str(i)],p.m);
%     V = subs(V,['l_' num2str(i)],p.l);
% end

% eq = solve(V==0,[q_1,q_2]);
%%
r = 2;
steps = 300;
theta = linspace(0,2*pi,steps);
ic = [r*cos(theta); r*sin(theta); zeros(1,steps); zeros(1,steps); zeros(1,steps); zeros(1,steps)]';
p.tspan = linspace(0,-30,1000);
for i = 1:steps
    p.q0 = ic(i,:);
    [t_tmp,q_tmp] = solveODE(p);
    Q(:,:,i) = q_tmp;
end
% hold on;
% for i = 1:steps
% plot(Q(:,1,i),Q(:,2,i))
% end
%%
hold on;
for i = 1:steps
% plot3(mod(Q(:,1,i),2*pi),mod(Q(:,2,i),2*pi),mod(Q(:,3,i),2*pi))
m1 = mod(Q(:,1,i),2*pi);
m2 = mod(Q(:,2,i),2*pi);
m3 = mod(Q(:,3,i),2*pi);
eps = 1;
m1(m1<eps)=NaN;
m2(m2<eps)=NaN;
m3(m3<eps)=NaN;
patch(m1,m2,m3,mod(Q(:,3,i),2*pi),'FaceColor','none','edgealpha',0.1)
drawnow;
pause(0.05)
end

%%
hold on;
for i = 1:steps
plot3(Q(:,2,i),sin(Q(:,1,i)),-cos(Q(:,1,i)))
drawnow;
view(50,30)
pause(0.05)
end