function u = controller(t,q,p,varargin)
ip=inputParser;
addOptional(ip,'traj',[p.q_des p.qd_des]);
addOptional(ip,'tspan',p.tspan);
parse(ip,varargin{:});
traj = ip.Results.traj;
tspan = ip.Results.tspan;

% int = interp1(tspan,traj,t);

% u = PD(q',p,int(1:p.n),int((1+p.n):2*p.n));
% 
% u = u.* p.control_action;
% u = Barrier(q,p);
u=zeros(1,p.n)';
end

function u = PD(q,p,q_des,qd_des)
u = p.kp.*(q(:,1:p.n)-q_des) + p.kd.*(q(:,(p.n+1):2*p.n)-qd_des);
end

function u = Barrier(q,p)
u = 1./q(1:p.n);
end

%     for i = 1:n
%         fprintf(fh,['res(' num2str(i+p.n) ') =res(' num2str(i+p.n) ')'...
%             ' + p.p_gain(' num2str(i) ')*(q( ' num2str(i) ')-p.q_des( ' num2str(i) ')) + p.d_gain(' num2str(i) ')*(q( ' num2str(i+p.n) ')) ;\n']);
%     end