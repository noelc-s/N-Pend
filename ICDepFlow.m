p = params();
p.redoEOM = false;

DeriveEOM(p);
[t,q] = solveODE(p);

clear state;

state(:,1,:) = q(:,1:3);
%% IC
res = 10;
[x,y,z] = sphere(res);
epsilon = 0.01;
x = x*epsilon;
y = y*epsilon;
z = z*epsilon;

for i = 1:numel(x)
    i
    q0 = [p.q0(1:3)+[x(i) y(i) z(i)]'; p.q0(4:6)];
    
    n = p.n;
    tspan = p.tspan;
    
    dq0 = zeros(2*n,1);
    
    [y0_new,yp0_new] = decic(@(t,q,dq) p.pendODE_i(t,q,dq),0,q0,[],dq0,[]);
    
    options = odeset('RelTol',.1,'AbsTol',.1);
    
    [t,q] = ode15i(p.pendODE_i,tspan,y0_new,yp0_new,options);
    state(:,i+1,:) = q(:,1:3);
end

csvwrite('t1.csv',state(:,:,1))
csvwrite('t2.csv',state(:,:,2))
csvwrite('t3.csv',state(:,:,3))
%%
% hold on;
% for i = 1:numel(x)
%     plot3(state(:,i,1),state(:,i,2),state(:,i,3));
% end

%%

x_axis_u = max(max(state(:,:,1)));
y_axis_u = max(max(state(:,:,2)));
z_axis_u = max(max(state(:,:,3)));
x_axis_l = min(min(state(:,:,1)));
y_axis_l = min(min(state(:,:,2)));
z_axis_l = min(min(state(:,:,3)));

figure('units','normalized','outerposition',[0 0 1 1])

for i = length(t):length(t)
    x_ = state(1:i,:,1);
    y_ = state(1:i,:,2);
    z_ = state(1:i,:,3);
    
    x_tmp=x_(end,:);
    y_tmp=y_(end,:);
    z_tmp=z_(end,:);
    
    x_(end,:)=NaN;
    y_(end,:)=NaN;
    z_(end,:)=NaN;
    
    scatter3(x_tmp',y_tmp',z_tmp',10,'filled','b');
    hold on;
    
    c = repmat((linspace(0,i/length(t),size(x_,1)))',1,size(x_,2));
    colormap(jet)
    patch(x_(:),y_(:),z_(:),c(:),'FaceColor','none','EdgeColor','interp')
    
    %     plot3(x_(:),y_(:),z_(:),'b');
    
    hold off;
    axis([x_axis_l,x_axis_u,...
        y_axis_l,y_axis_u,...
        z_axis_l,z_axis_u]);
    view(i/length(t)*50,30)

    
    drawnow
    M(i) = getframe(gcf);
    pause(0.01);
end
%%
video = VideoWriter('tmp.avi','Uncompressed AVI');
video.FrameRate = 60;
open(video);
writeVideo(video,M);
close(video);