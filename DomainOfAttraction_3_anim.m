clear;close all;clc;
addpath(genpath(pwd));
p = params();
DeriveEOM(p);


[~,E] = DeriveEnergy(p);

[X,Y,Z]=meshgrid(-2*pi:0.5:2*pi,-2*pi:0.5:2*pi,-2*pi:0.5:6*pi);

[t,q] = solveODE(p);
%%
for k = 1:length(t)
    e(k) = p.pendEnergy(q(k,:),p);
end

figure('units','normalized','outerposition',[0 0 1 1])
% figure
for g = 1:length(t)
    
    for i = 1:size(X,1)
        for j = 1:size(Y,2)
            for k = 1:size(Z,3)
                en(i,j,k)=p.pendEnergy([X(i,j,k),Y(i,j,k),Z(i,j,k),q(g,4),q(g,5),q(g,6)],p);
            end
        end
    end
    
    hold on
    for i = 1:size(X,3)
        x = X(:,:,i);
        y = Y(:,:,i);
        z = Z(:,:,i);
        e = en(:,:,i);
        scatter3(x(:),y(:),z(:),200./e(:),'b','filled')
    end
    plot3(q(1:g,1),q(1:g,2),q(1:g,3),'LineWidth',3)
    scatter3(q(g,1),q(g,2),q(g,3),100,'filled','LineWidth',3)
    view(g/10+30,30)
    axis off;
    drawnow;
    M(g) = getframe(gcf);
    clf;
%     pause
%     
end

%% 
video = VideoWriter('tmp.avi','Motion JPEG AVI');
video.FrameRate = 24;
open(video);
writeVideo(video,M);
close(video);


% %
% %
% 
% % for i = 1:length(t)
% %     e(i) = p.pendEnergy(q(i,:),p);
% % end
% %
% hold on;
% % surf(X,Y,Z);
% plot3(q(:,1),q(:,2),q(:,3),'LineWidth',2)
% % contour(X,Y,Z,[1:25])