clear;close all;clc;
addpath(genpath(pwd));
p = params();
sol = DeriveEOM(p);
tdd = sym('A',[1 p.n]);
%%
for i = 1:p.n
    tdd(i) = eval(['sol.tdd_' num2str(i)']);
end

vars = {'m','t','td','tdd','l','b','u'};

for i = 1:p.n
    for var=vars
        eval(['syms ' var{:} '_' num2str(i)]);
        assume(eval([var{:} '_' num2str(i)]),'real');
    end
end

for i = 1:p.n
    for j = 1:p.n
        tdd(i) = subs(tdd(i),eval(['u_' num2str(j)]),0);
        tdd(i) = subs(tdd(i),eval(['b_' num2str(j)]),p.b(j));
        tdd(i) = subs(tdd(i),eval(['l_' num2str(j)]),p.l(j));
        tdd(i) = subs(tdd(i),eval(['m_' num2str(j)]),p.m(j));
    end
end

state_vars = {'t','td'};

q_var = [];

for var=state_vars
    for i = 1:p.n
        q_var = [q_var eval([var{:} '_' num2str(i)])];
    end
end


tdd(p.n+1:2*p.n) = tdd(1:p.n);
tdd(1:p.n) = q_var(p.n+1:2*p.n)';

for i = 1:2*p.n
    for j = 1:2*p.n
        Dxf(i,j) = diff(tdd(i),q_var(j));
    end
end
Dxf0 = eval(subs(Dxf,{'t_1','t_2','t_3','td_1','td_2','td_3'},{0,0,0,0,0,0}));
eig(Dxf0)

Dxf_tau = repmat(tdd',1,2*p.n);

for i = 1:p.n
    for j = 1:2*p.n
        for k = 1:2*p.n
            if k~=j
                Dxf_tau(i,j) = subs(Dxf_tau(i,j),q_var(k),0);
            end
        end
        Dxf_tau(i,j) = subs(Dxf_tau(i,j),q_var(j),1);
    end
end

for i = p.n+1:2*p.n
    for j = 1:2*p.n
        for k = 1:2*p.n
            if k~=j
                Dxf_tau(i,j) = subs(Dxf_tau(i,j),q_var(k),0);
            end
        end
        Dxf_tau(i,j) = Dxf_tau(i,j)/q_var(j);
    end
end

%%
P = lyap(Dxf0, eye(2*p.n));
lambda_P = eig(P);

N_max = 1/(2*max(lambda_P))

G = Dxf_tau - Dxf0;
G_ = matlabFunction(G);
vars = regexp(cellstr(sprintf('a(%d) ',1:nargin(G_))),' ','split');
vars = strjoin(vars{1},', ');
G_norm = eval(['@(a) -norm(G_(' vars '))']);
options = optimoptions('fmincon','Display','iter','StepTolerance',1e-20);
[x,fval] = fmincon(G_norm,0.01*ones(1,nargin(G_)),[],[],[],[],-100*ones(1,nargin(G_)),100*ones(1,nargin(G_)),@(x)nonlcon(x,G_norm,N_max),options)


function [c,ceq] = nonlcon(x,G_norm,N_max)
ceq = 0;
c = -G_norm(x) - N_max;
end

