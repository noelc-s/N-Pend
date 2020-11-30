%% Explicit Form
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
    tdd = subs(tdd,['b_' num2str(i)],b);
    tdd = subs(tdd,['m_' num2str(i)],m);
    tdd = subs(tdd,['l_' num2str(i)],l);
end

J = [zeros(length(tdd),length(tdd)) eye(length(tdd));jacobian(tdd)];
J_ = matlabFunction(J);
vars = regexp(cellstr(sprintf('a(%d) ',1:nargin(J_))),' ','split');
vars = strjoin(vars{1},', ');
J_norm = eval(['@(a) -norm(J_(' vars '))']);

[x,fval] = fmincon(J_norm,ones(1,nargin(J_)),[],[],[],[],-1*ones(1,nargin(J_)),1*ones(1,nargin(J_)))
% [x,fval] = fmincon(J_norm,ones(1,nargin(J_)),[],[],[],[])

%%
L_min = inf;
L_max = -inf;

for i = 1:size(q,2)
    for j = 1:size(q,1)
        for k = j:size(q,1)
            [L_min,tmp] = min([L_min,(q(j,i)-q(k,i))/(t(j)-t(k))]);
            if tmp == 2
                min_ind = [i j];
            end
            [L_max,tmp] = max([L_max,(q(j,i)-q(k,i))/(t(j)-t(k))]);
            if tmp == 2
                max_ind = [i j];
            end 
        end
    end
end
%%
plot(t,q)
hold on;

x = t;
y_u = L_max*(t-t(max_ind(2)))+q(max_ind(2),max_ind(1));
y_d = -L_max*(t-t(max_ind(2)))+q(max_ind(2),max_ind(1));
plot(x,y_u)
plot(x,y_d)