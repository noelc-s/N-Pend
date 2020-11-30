function dq = pendOde(t,q,p)
t
dq = q;
t_1=q(1,:);
td_1=q(3,:);
t_2=q(2,:);
td_2=q(4,:);
u=p.controller(t,q,p);
u_1=u(1,:);
u_2=u(2,:);
dq(1,:)=q(3,:);
dq(2,:)=q(4,:);
dq(3,:)=-(100.*p.l(1).*u_2.*cos(t_1 - t_2) - 100.*p.l(2).*u_1 - 100.*p.b(1).*p.l(2).*td_1 + 981.*p.l(1).*p.l(2).*p.m(1).*sin(t_1) + 981.*p.l(1).*p.l(2).*p.m(2).*sin(t_1) + 100.*p.b(2).*p.l(1).*td_2.*cos(t_1 - t_2) - 981.*p.l(1).*p.l(2).*p.m(2).*cos(t_1 - t_2).*sin(t_2) + 100.*p.l(1).*p.l(2).^2.*p.m(2).*td_2.^2.*sin(t_1 - t_2) + 100.*p.l(1).^2.*p.l(2).*p.m(2).*td_1.^2.*cos(t_1 - t_2).*sin(t_1 - t_2))/(100.*p.l(1).^2.*p.l(2).*(p.m(1) + p.m(2) - p.m(2).*cos(t_1 - t_2).^2));
dq(4,:)=(100.*p.l(1).*p.m(1).*u_2 + 100.*p.l(1).*p.m(2).*u_2 + 100.*p.b(2).*p.l(1).*p.m(1).*td_2 + 100.*p.b(2).*p.l(1).*p.m(2).*td_2 - 100.*p.l(2).*p.m(2).*u_1.*cos(t_1 - t_2) - 981.*p.l(1).*p.l(2).*p.m(2).^2.*sin(t_2) - 100.*p.b(1).*p.l(2).*p.m(2).*td_1.*cos(t_1 - t_2) + 100.*p.l(1).^2.*p.l(2).*p.m(2).^2.*td_1.^2.*sin(t_1 - t_2) + 981.*p.l(1).*p.l(2).*p.m(2).^2.*cos(t_1 - t_2).*sin(t_1) - 981.*p.l(1).*p.l(2).*p.m(1).*p.m(2).*sin(t_2) + 981.*p.l(1).*p.l(2).*p.m(1).*p.m(2).*cos(t_1 - t_2).*sin(t_1) + 100.*p.l(1).*p.l(2).^2.*p.m(2).^2.*td_2.^2.*cos(t_1 - t_2).*sin(t_1 - t_2) + 100.*p.l(1).^2.*p.l(2).*p.m(1).*p.m(2).*td_1.^2.*sin(t_1 - t_2))/(100.*p.l(1).*p.l(2).^2.*p.m(2).*(p.m(1) + p.m(2) - p.m(2).*cos(t_1 - t_2).^2));
dq(4,:)=dq(4,:)+controller(t,p,q);
end
function u = controller(t,p,q)
t_1=q(1,:);
td_1=q(3,:);
t_2=q(2,:);
td_2=q(4,:);

fk = [sin(t_1)+sin(t_2) -cos(t_1)-cos(t_2)];
% h = 1-fk(1)
h = 1-sin(t_1)-sin(t_2);
Jh = [-cos(t_1) -cos(t_2)];

alpha = 10;
alpha_e = 10;

D = [(p.m(1)+p.m(2))*p.l(1)^2+p.m(2)*p.l(2)^2+2*p.m(2)*p.l(1)*p.l(2)*cos(t_2-t_1) p.m(2)*p.l(2)^2+p.m(2)*p.l(1)*p.l(2)*cos(t_2-t_1);p.m(2)*p.l(2)^2+p.m(2)*p.l(1)*p.l(2)*cos(t_2-t_1) p.m(2)*p.l(2)^2];
C = [0 -p.m(2)*p.l(1)*p.l(2)*(2*td_1+(td_2-td_1))*sin(td_2-td_1); p.m(2)*p.l(1)*p.l(2)*td_1*sin(t_2-t_1) 0];
G = p.g*[(p.m(1)+p.m(2))*p.l(1)*sin(t_1)+p.m(2)*p.l(2)*sin(t_2); p.m(2)*p.l(2)*sin(t_2)];
    % Is this right?
    % https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-832-underactuated-robotics-spring-2009/readings/MIT6_832s09_read_appA.pdf
h_D = -0.5*[td_1 td_2]*D*[td_1; td_2]+alpha_e*h;
B = [0;1];

u_des=sin(10*t);

Q = 1;
R = -2*u_des;
A = [td_1 td_2]*B;
b = alpha*h_D + alpha_e*Jh*[td_1; td_2] + G'*[td_1; td_2];

u=quadprog(Q,R,A,b)
% u=u_des;

end