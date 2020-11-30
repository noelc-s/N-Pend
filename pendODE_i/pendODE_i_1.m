function res = pendOde(t,q,dq,p)
res = zeros(2,1);
res(1) =dq(1) - q(2);
res(2) =(981*p.l(1)*p.m(1)*sin(q(1)))/100 - p.b(1)*q(2) + dq(2)*p.l(1)^2*p.m(1);
u=p.controller(t,q,p);
res(2) =res(2) + u(1);
end