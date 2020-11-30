clear;close all;clc;
addpath(genpath(pwd));
p = params();
p.redoEOM = true;
p.explicit = true;
p.controls = true;

f = DeriveEOM(p);


[~,V,dV] = DeriveEnergy(p,'equationForm',true);

vars = {'q','tdd','b','m','l'};

for i = 1:2*p.n
    for var=vars
        eval(['syms ' var{:} '_' num2str(i)]);
        eval(['assumeAlso(' var{:} '_' num2str(i) ', ''real'');']);

    end
end

for i = 1:p.n
    V = subs(V,['b_' num2str(i)],p.b);
    V = subs(V,['m_' num2str(i)],p.m);
    V = subs(V,['l_' num2str(i)],p.l);
end

%%
q_dd = [];
for i = 1:p.n
    eval(['q_dd = [q_dd q_' num2str(i)+p.n '];']); 
end
fn = fields(f);
for j = 1:length(fn)
    for i = 1:p.n
        eval(['f.' fn{j} '= subs(f.' fn{j} ',[''t_'' num2str(i)],[''q_'' num2str(i)]);']);
        eval(['f.' fn{j} '= subs(f.' fn{j} ',[''td_'' num2str(i)],[''q_'' num2str(i)+p.n]);']);
    end
    q_dd = [q_dd f.(fn{j})];
    eval(['assumeAlso(' fn{j} ', ''real'');']);
end

V_dot = simplify(dV*q_dd')
