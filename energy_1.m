function [V,Vd] = energy(q,p)
V =(p.l(1)^2*p.m(1)*q(2)^2)/2 - (981*p.l(1)*p.m(1)*cos(q(1)))/100;
Vd(1) =(981*p.l(1)*p.m(1)*sin(q(1)))/100;
end