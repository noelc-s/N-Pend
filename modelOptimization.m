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
x_0 = XE_MLF(p.q0(1),p.q0(2));
y_0 = YE_MLF(p.q0(1),p.q0(2));

ord = 2;

x=getBezier(ord,[x_0 p.n*rand(1,ord)-p.n/2]);
y=getBezier(ord,[y_0 p.n*rand(1,ord)-p.n/2]);
% x = @(t) sum(p.l)*cos(t);
% y = @(t) -sum(p.l)*sin(t);
% x=@(t) t-1.5;
% y = @(t) 0*t-1;

tau = linspace(0,1);

%%
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','StepTolerance',1e-10,'FunctionTolerance',1e-10,'RelLineSrchBnd',1e-10,'UseParallel',true);
[par,val] = fmincon(@(par) getMetric(par,p,x,y,tau,XE_MLF,YE_MLF),2*ones(1,p.n),[],[],[],[],-30*ones(1,p.n),30*ones(1,p.n),[],options)
%%
p.q0(p.n+1:2*p.n) = par(1:p.n)';
[t,q]=solveODE(p);
plot(x(tau),y(tau))
animatePend
%%

function metric = getMetric(par,p,x,y,tau,XE_MLF,YE_MLF)
p.q0(p.n+1:2*p.n) = par(1:p.n)';

[t,q]=solveODE(p);

X = XE_MLF(q(:,1),q(:,2));
Y = YE_MLF(q(:,1),q(:,2));

x_err = norm(X-x(tau)',1);
y_err = norm(Y-y(tau)',1);
metric = norm([x_err y_err],1);
end
%%
function B = getBezier(n,coefs)
B = @(t) (1-t).^n*coefs(1);
for i = 1:n
    B = @(t) B(t) + nchoosek(n,i).*(1-t).^(n-i).*t.^i*coefs(i+1);
end
end