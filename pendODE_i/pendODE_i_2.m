function res = pendOde(t,q,dq,p)
res = zeros(4,1);
res(1) =dq(2) - q(4);
res(2) =dq(1) - q(3);
res(3) =(981*p.l(1)*p.m(1)*sin(q(1)))/100 - p.b(1)*q(3) + (981*p.l(1)*p.m(2)*sin(q(1)))/100 + dq(3)*p.l(1)^2*p.m(1) + dq(3)*p.l(1)^2*p.m(2) + p.l(1)*p.l(2)*p.m(2)*(dq(4)*cos(q(1) - q(2)) - q(4)*sin(q(1) - q(2))*(q(3) - q(4))) + p.l(1)*p.l(2)*p.m(2)*q(3)*q(4)*sin(q(1) - q(2));
res(4) =(981*p.l(2)*p.m(2)*sin(q(2)))/100 - p.b(2)*q(4) + dq(4)*p.l(2)^2*p.m(2) + p.l(1)*p.l(2)*p.m(2)*(dq(3)*cos(q(1) - q(2)) - q(3)*sin(q(1) - q(2))*(q(3) - q(4))) - p.l(1)*p.l(2)*p.m(2)*q(3)*q(4)*sin(q(1) - q(2));
u=p.controller(t,q,p);
res(3) =res(3) + u(1);
res(4) =res(4) + u(2);
end