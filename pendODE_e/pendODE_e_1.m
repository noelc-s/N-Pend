function dq = pendOde(t,q,p)
dq = q;
t_1=q(1,:);
td_1=q(2,:);
u=p.controller(t,q,p);
u_1=u(1,:);
dq(1,:)=q(2,:);
dq(2,:)=(u_1 + p.b(1).*td_1 - (981.*p.l(1).*p.m(1).*sin(t_1))/100)/(p.l(1).^2.*p.m(1));
end