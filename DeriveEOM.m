function [sol] = DeriveEOM(p)
if p.redoEOM | ~exist(['pendODE_i/pendODE_i_' num2str(p.n) '.m']) | (~exist(['pendODE_e/pendODE_e_' num2str(p.n) '.m']) & p.explicit)
    disp('Deriving Equations of Motion')
    n = p.n;
    g = p.g;
    m = p.m;
    b = p.b;
    l = p.l;
    
    vars = {'m','t','td','tdd','l','b','u'};
    
    for i = 1:n
        for var=vars
            eval(['syms ' var{:} '_' num2str(i)]);
        end
    end
    
    T = 0;
    V = 0;
    L = 0;
    
    
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
    disp(' Lagrangian Formed')
    
    %%
    
    delta = @(a,b) a==b;
    
    dL_dt = sym('A',[1 n]);
    ddL_dtd_dt = sym('A',[1 n]);
    
    for i = 1:n
        dL_dt(i) = 0;
        for j = i:n
            for var=vars
                eval([var{1} '_i = eval([''' var{1} '_'' num2str(i)]);']);
                eval([var{1} '_j = eval([''' var{1} '_'' num2str(j)]);']);
                eval([var{1} '_k = eval([''' var{1} '_'' num2str(k)]);']);
            end
            dL_dt(i) = dL_dt(i) - g*m_j*l_i*sin(t_i);
        end
    end
    for i = 1:n
        for j = 1:n
            for k = max(i,j):n
                for var=vars
                    eval([var{1} '_i = eval([''' var{1} '_'' num2str(i)]);']);
                    eval([var{1} '_j = eval([''' var{1} '_'' num2str(j)]);']);
                    eval([var{1} '_k = eval([''' var{1} '_'' num2str(k)]);']);
                end
                dL_dt(i) = dL_dt(i) + l_i*td_i*l_j*td_j*m_k*sin(t_j-t_i)*(1-delta(i,j));
            end
        end
    end
    
    for i = 1:n
        ddL_dtd_dt(i) = 0;
        for j = i:n
            for var=vars
                eval([var{1} '_i = eval([''' var{1} '_'' num2str(i)]);']);
                eval([var{1} '_j = eval([''' var{1} '_'' num2str(j)]);']);
                eval([var{1} '_k = eval([''' var{1} '_'' num2str(k)]);']);
            end
            ddL_dtd_dt(i) = ddL_dtd_dt(i)+ l_i^2*tdd_i*m_j;
        end
    end
    
    for i = 1:n
        for j = 1:n
            for k = max(i,j):n
                for var=vars
                    eval([var{1} '_i = eval([''' var{1} '_'' num2str(i)]);']);
                    eval([var{1} '_j = eval([''' var{1} '_'' num2str(j)]);']);
                    eval([var{1} '_k = eval([''' var{1} '_'' num2str(k)]);']);
                end
                ddL_dtd_dt(i) = ddL_dtd_dt(i) + l_i*l_j*m_k*(tdd_j*cos(t_i-t_j)-td_j*(td_i-td_j)*sin(t_i-t_j))*(1-delta(i,j));
            end
        end
    end
    
    EL = sym('A',[1 n]);
    EL_imp = sym('A',[1 n]);
    
    for i = 1:n
        for var=vars
            eval([var{1} '_i = eval([''' var{1} '_'' num2str(i)]);']);
            eval([var{1} '_j = eval([''' var{1} '_'' num2str(j)]);']);
            eval([var{1} '_k = eval([''' var{1} '_'' num2str(k)]);']);
        end
        EL(i) = ddL_dtd_dt(i) - dL_dt(i) == b_i*td_i + u_i;
        EL_imp(i) = ddL_dtd_dt(i) - dL_dt(i) - b_i*td_i;
    end
    
    disp(' Euler-Lagrange Equations Written')
    
    %%
    
    Q = {'q','dq'};
    
    for i = 1:2*n
        for q=Q
            eval(['syms ' q{:} '_' num2str(i)]);
        end
    end
    
    for i = 1:n
        EL_imp = subs(EL_imp,['t_' num2str(i)],['q_' num2str(i)]);
        EL_imp = subs(EL_imp,['td_' num2str(i)],['q_' num2str(i+n)]);
        EL_imp = subs(EL_imp,['tdd_' num2str(i)],['dq_' num2str(i+n)]);
    end
    
    for i = 1:n
        EL_imp = eval(['[dq_' num2str(i) '-q_' num2str(i+n) ' EL_imp]']);
    end
    disp(' Euler-Lagrange Equations Solved')
    %%
    if isempty(p.filename_i)
        fh = fopen(['pendODE_i/pendODE_i_' num2str(n) '.m'],'w');
    else
        fh = fopen(['pendODE_i/' p.filename_i num2str(n) '.m'],'w');
    end
    fprintf(fh,'function res = pendOde(t,q,dq,p)\n');
    fprintf(fh,['res = zeros(' num2str(2*n) ',1);\n']);
    for i = 1:2*n
        eq = char(EL_imp(i));
        for j = 2*n:-1:1
            eq = replace(eq,['b_' num2str(j)],['p.b(' num2str(j) ')']);
            eq = replace(eq,['l_' num2str(j)],['p.l(' num2str(j) ')']);
            eq = replace(eq,['m_' num2str(j)],['p.m(' num2str(j) ')']);
            eq = replace(eq,['q_' num2str(j)],['q(' num2str(j) ')']);
            eq = replace(eq,['dq_' num2str(j)],['dq(' num2str(j) ')']);
        end
        fprintf(fh,['res(' num2str(i) ') =' eq ';\n']);
    end
    if p.controls
    fprintf(fh,'u=p.controller(t,q,p);\n');
    for i = 1:n
        fprintf(fh,['res(' num2str(i+p.n) ') =res(' num2str(i+p.n) ') + u(' num2str(i) ');\n']);
    end
    end
    fprintf(fh,'end');
    fclose(fh);
    disp(' Implicit Form Written')
    %%
    if p.explicit
        tdd = [];
        for i = 1:n
            tdd = [tdd eval(['tdd_' num2str(i)])];
        end
        
        sol = solve(EL,tdd);
        
        if ~isstruct(sol)
            tmp = sol; clear sol; sol.tdd_1 = tmp;
        end
        
        tdd = sym('A',[1 n]);
        for i = 1:n
            tdd(i) = eval(['sol.tdd_' num2str(i)']);
        end
        
        for i = 1:n
            c_tdd{i} = char(tdd(i));
            c_tdd{i} = replace(c_tdd{i},'*' ,'.*');
            c_tdd{i} = replace(c_tdd{i},'^' ,'.^');
            for j = 1:n
            c_tdd{i} = replace(c_tdd{i},['b_' num2str(j)],['p.b(' num2str(j) ')']);
            c_tdd{i} = replace(c_tdd{i},['l_' num2str(j)],['p.l(' num2str(j) ')']);
            c_tdd{i} = replace(c_tdd{i},['m_' num2str(j)],['p.m(' num2str(j) ')']);
            end
        end
        
        if isempty(p.filename_e)
            fh = fopen(['pendODE_e/pendODE_e_' num2str(n) '.m'],'w');
        else
            fh = fopen(['pendODE_e/' p.filename_e num2str(n) '.m'],'w');
        end
        fprintf(fh,'function dq = pendOde(t,q,p)\n');
        fprintf(fh,'dq = q;\n');
        for i = 1:n
            fprintf(fh,['t_' num2str(i) '=q(' num2str(i) ',:);\n']);
            fprintf(fh,['td_' num2str(i) '=q(' num2str(i+n) ',:);\n']);
        end
        fprintf(fh,'u=p.controller(t,q,p);\n');
        for i = 1:n
            fprintf(fh,['u_' num2str(i) '=u(' num2str(i) ',:);\n']);
        end
        for i = 1:n
            fprintf(fh,['dq(' num2str(i) ',:)=q(' num2str(i+n) ',:);\n']);
        end
        for i = 1:n
            fprintf(fh,['dq(' num2str(i+n) ',:)=' c_tdd{i} ';\n']);
        end
        fprintf(fh,'end');
        fclose(fh);
        disp(' Explicit Form Written')
    end
end

end
%% APPENDIX
%% ode45
% This is an equivalent way to solve for equations of motion, but solve was
% very slow
%
%
% tdd = [];
% for i = 1:n
%     tdd = [tdd eval(['tdd_' num2str(i)])];
% end
%
% sol = solve(EL,tdd);
%
% if ~isstruct(sol)
%     tmp = sol; clear sol; sol.tdd_1 = tmp;
% end
%
% tdd = sym('A',[1 n]);
% for i = 1:n
%     tdd(i) = eval(['sol.tdd_' num2str(i)']);
% end
%
% for i = 1:n
%     tdd = subs(tdd,['b_' num2str(i)],b);
%     tdd = subs(tdd,['m_' num2str(i)],m);
%     tdd = subs(tdd,['l_' num2str(i)],l);
% end
% disp('Section 3 Complete')
%
% fh = fopen('pendODE.m','w');
% fprintf(fh,'function dq = pendOde(t,q)\n');
% fprintf(fh,'dq = q;\n');
% for i = 1:n
%     fprintf(fh,['t_' num2str(i) '=q(' num2str(i) ');\n']);
%     fprintf(fh,['td_' num2str(i) '=q(' num2str(i+n) ');\n']);
% end
% for i = 1:n
%     fprintf(fh,['dq(' num2str(i) ')=q(' num2str(i+n) ');\n']);
% end
% for i = 1:n
%     fprintf(fh,['dq(' num2str(i+n) ')=' char(tdd(i)) ';\n']);
% end
% fprintf(fh,'end');
% fclose(fh);
%
% tspan = linspace(0,10,1000);
% q0 = [ones(1,n) zeros(1,n)];
%
% [t,q] = ode45(@(t,q) pendODE(t,q),tspan,q0);
%
% %%
% y1 = linspace(-2,8,100);
% y2 = linspace(-8,8,100);
% [x,y] = meshgrid(y1,y2);
% u = zeros(size(x));
% v = zeros(size(x));
% t=0; % we want the derivatives at each point at t=0, i.e. the starting time
%
% link = 1;
%
% for i = 1:numel(x)
%     ic = [];
%     for j = 1:n
%         if j == link
%             ic = [ic; x(i)];
%         else
%             ic = [ic; 0];
%         end
%     end
%     for j = 1:n
%         if j == link
%             ic = [ic; y(i)];
%         else
%             ic = [ic; 0];
%         end
%     end
%     Yprime = pendODE(t,ic);
%     u(i) = Yprime(link);
%     v(i) = Yprime(link+n);
% end
%
% % quiver3(cos(x),sin(x),y,cos(u),sin(u),v,'r'); figure(gcf)
% quiver(x,y,u,v,'r'); figure(gcf)
%
% hold on
% for y20 = 0:1:8
%     ic = zeros(n,1);
%     for j = 1:n
%         if j==link
%             ic = [ic; y20];
%         else
%             ic = [ic; 0];
%         end
%     end
%     [ts,ys] = ode45(@(t,q) pendODE(t,q),[0,3],ic);
%     plot(ys(:,link),ys(:,link+n))
%     plot(ys(1,link),ys(1,link+n),'bo') % starting point
%     plot(ys(end,link),ys(end,link+n),'ks') % ending point
% end
% hold off
