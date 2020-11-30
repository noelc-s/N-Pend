clear;close all;clc;
addpath(genpath(pwd));
p = params();
DeriveEOM(p);


[~,E] = DeriveEnergy(p);

[X,Y]=meshgrid(-2*pi:0.1:2*pi);

[t,q] = solveODE(p);

for k = 1:length(t)
    e(k) = p.pendEnergy(q(k,:),p);
end
% 
% plot3(q(:,1),q(:,2),e(:))

%%

figure('units','normalized','outerposition',[0 0 1 1])
for k = 1:length(t)
%     e(k) = p.pendEnergy(q(k,:),p);
    
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            Z(i,j)=p.pendEnergy([X(i,j),Y(i,j),q(k,3),q(k,4)],p);
        end
    end
    hold on;
    surf(X,Y,Z,'FaceAlpha',0.2,'EdgeAlpha',0.1);
    plot3(q(1:k,1),q(1:k,2),e(1:k),'LineWidth',3)
    scatter3(q(k,1),q(k,2),e(k),100,'filled','LineWidth',3)
%     contour(X,Y,Z,[1:100])
%     plot(q(1:k,1),q(1:k,2),'LineWidth',2);
    view(10,60)
    axis([-2*pi 2*pi -2*pi 2*pi 0 100])
    axis off
    hold off;
    drawnow
        M(k) = getframe(gcf);

    clf;
    
end

%%
% for i = 1:2:length(M)
%    M_tmp(i/2+0.5) = M(i);
% end

video = VideoWriter('tmp.avi','Motion JPEG AVI');
video.FrameRate = 24;
open(video);
writeVideo(video,M);
close(video);

