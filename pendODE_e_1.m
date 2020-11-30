function dq = pendOde(t,q,p)
dq = q;
t_1=q(1,:);
td_1=q(2,:);
u=p.controller(t,q,p);
u_1=u(1,:);
dq(1,:)=q(2,:);
dq(2,:)=u_1 - td_1/10 - (981.*sin(t_1))/100;
end