function [L,E,dE] = DeriveEnergy(p,varargin)
ip=inputParser;
addOptional(ip,'equationForm','false',@(x)islogical(x));
parse(ip,varargin{:});
eF = ip.Results.equationForm;

n = p.n;
g = p.g;
m = p.m;
b = p.b;
l = p.l;

vars = {'m','t','td','tdd','l','b'};

for i = 1:n
    for var=vars
        eval(['syms ' var{:} '_' num2str(i)]);
    end
end

T = 0;
V = 0;
L = 0;
E = 0;


for i = 1:n
    for j = 1:i
        for k = 1:i
            for var=vars
                eval([var{1} '_i = eval([''' var{1} '_'' num2str(i)]);']);
                eval([var{1} '_j = eval([''' var{1} '_'' num2str(j)]);']);
                eval([var{1} '_k = eval([''' var{1} '_'' num2str(k)]);']);
            end
            T = T + 1/2*m_i*l_j*l_k*td_j*td_k*cos(t_j-t_k);
        end
    end
end
for i = 1:n
    for j = 1:i
        for var=vars
            eval([var{1} '_i = eval([''' var{1} '_'' num2str(i)]);']);
            eval([var{1} '_j = eval([''' var{1} '_'' num2str(j)]);']);
            eval([var{1} '_k = eval([''' var{1} '_'' num2str(k)]);']);
        end
        V = V - g*m_i*l_j*cos(t_j);
    end
end

L = L+(T-V);
E = E+T+V;

for i = 1:p.n
    dE(i) = eval(['diff(E, t_' num2str(i) ')']);
    dE(i+p.n) = eval(['diff(E, td_' num2str(i) ')']);
end

%%

if eF
    if p.redoEnergy | ~exist(['Energy/energy_' num2str(p.n) '.m'])
        for i = 1:n
            E = subs(E,['t_' num2str(i)],['q_' num2str(i)]);
            E = subs(E,['td_' num2str(i)],['q_' num2str(i+n)]);
            E = subs(E,['tdd_' num2str(i)],['dq_' num2str(i+n)]);
            dE = subs(dE,['t_' num2str(i)],['q_' num2str(i)]);
            dE = subs(dE,['td_' num2str(i)],['q_' num2str(i+n)]);
            dE = subs(dE,['tdd_' num2str(i)],['dq_' num2str(i+n)]);
        end
        if isempty(p.energyName)
            fh = fopen(['Energy/energy_' num2str(n) '.m'],'w');
        else
            fh = fopen(['Energy/' p.energyName num2str(n) '.m'],'w');
        end
        fprintf(fh,'function [V,Vd] = energy(q,p)\n');
        
        eq = char(E);
        for j = 2*n:-1:1
            eq = replace(eq,['b_' num2str(j)],['p.b(' num2str(j) ')']);
            eq = replace(eq,['l_' num2str(j)],['p.l(' num2str(j) ')']);
            eq = replace(eq,['m_' num2str(j)],['p.m(' num2str(j) ')']);
            eq = replace(eq,['q_' num2str(j)],['q(' num2str(j) ')']);
            eq = replace(eq,['dq_' num2str(j)],['dq(' num2str(j) ')']);
        end
        fprintf(fh,['V =' eq ';\n']);
        for i = 1:p.n
            eq = char(dE(i));
            for j = 2*n:-1:1
                eq = replace(eq,['b_' num2str(j)],['p.b(' num2str(j) ')']);
                eq = replace(eq,['l_' num2str(j)],['p.l(' num2str(j) ')']);
                eq = replace(eq,['m_' num2str(j)],['p.m(' num2str(j) ')']);
                eq = replace(eq,['q_' num2str(j)],['q(' num2str(j) ')']);
                eq = replace(eq,['dq_' num2str(j)],['dq(' num2str(j) ')']);
            end
            fprintf(fh,['Vd(' num2str(i) ') =' eq ';\n']);
        end
        fprintf(fh,'end');
        fclose(fh);
    end
end
end