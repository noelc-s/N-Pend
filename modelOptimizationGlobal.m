%%
p=params();
vars = {'th'};
th_arr = [];
for i = 1:p.n
    for var=vars
        eval(['syms ' var{:} '_' num2str(i)]);
        th_arr = [th_arr eval([var{:} '_' num2str(i)])];
    end
end

X_RHS = 0;
Y_RHS = 0;
for i = 1:p.n
X_RHS = X_RHS+p.l(i)*cos(eval(['th_' num2str(i)]));
Y_RHS = Y_RHS+p.l(i)*sin(eval(['th_' num2str(i)]));
end

XE_MLF = matlabFunction(Y_RHS,'Vars',th_arr);
YE_MLF = matlabFunction(-X_RHS,'Vars',th_arr);
%%
x_0 = XE_MLF(p.q0(1),p.q0(2),p.q0(3));
y_0 = YE_MLF(p.q0(1),p.q0(2),p.q0(3));

ord = 2;

% x=getBezier(ord,[x_0 p.n*rand(1,ord)-p.n/2]);
% y=getBezier(ord,[y_0 p.n*rand(1,ord)-p.n/2]);
% x = @(t) sum(p.l)*cos(t);
% y = @(t) -sum(p.l)*sin(t);
% x=@(t) t-1.5;
% y = @(t) 0*t-1;

x =@(t) -0.9093-0.5*sin(pi*t);
y =@(t) -1.0839-0.5*cos(pi*t);

tau = linspace(0,1);
% plot(x(tau),y(tau))
% axis equal

%%
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','StepTolerance',1e-10,'FunctionTolerance',1e-10,'RelLineSrchBnd',1e-10,'UseParallel',true);
% problem = createOptimProblem('fmincon','objective',@(par) getMetric(par,p,x,y,tau,XE_MLF,YE_MLF),'x0',2*ones(1,2*p.n),'options',options);
% [par,val] = fmincon(problem)
%%
% gs = GlobalSearch('Display','iter');
% rng(14,'twister');
% [par,fval] = run(gs,problem);
%% 
% options = optimoptions('fminunc','Display','iter','Algorithm','quasi-newton','StepTolerance',1e-10,'FunctionTolerance',1e-10,'UseParallel',true);
% problem = createOptimProblem('fminunc','objective',@(par) getMetric(par,p,x,y,tau,XE_MLF,YE_MLF),'x0',2*ones(1,p.n),'options',options);
% [par,val] = fminunc(problem)
%%
% ms = MultiStart;
% [par,val] = run(ms,problem,10)
%%
options = optimoptions('ga','Display','iter','UseParallel',true);
[par, val] = ga(@(par) getMetric(par,p,x,y,tau,XE_MLF,YE_MLF),6,[],[],[],[],[-10*ones(1,p.n) 0.5*ones(1,p.n)],[10*ones(1,p.n) 10*ones(1,p.n)],[],options);
%%
p.q0(p.n+1:2*p.n) = par(1:p.n)';
p.m = par(p.n+1:2*p.n);
[t,q]=solveODE(p);
plot(x(tau),y(tau))
animatePend
%%

function metric = getMetric(par,p,x,y,tau,XE_MLF,YE_MLF)
p.q0(p.n+1:2*p.n) = par(1:p.n)';
p.m = par(p.n+1:2*p.n);

[t,q]=solveODE(p);

X = XE_MLF(q(:,1),q(:,2),q(:,3));
Y = YE_MLF(q(:,1),q(:,2),q(:,3));

x_err = norm(X-x(tau)',2);
y_err = norm(Y-y(tau)',2);
metric = norm([x_err y_err],2);
end
%%
function B = getBezier(n,coefs)
B = @(t) (1-t).^n*coefs(1);
for i = 1:n
    B = @(t) B(t) + nchoosek(n,i).*(1-t).^(n-i).*t.^i*coefs(i+1);
end
end