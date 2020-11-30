clear;close all;clc;
addpath(genpath(pwd));
p = params();
DeriveEOM(p);
plotEnergy = false;


[~,V] = DeriveEnergy(p);
[t,q] = solveODE(p);



%% Plotting Pend

if plotEnergy
    [X,Y]=meshgrid(-2*pi:0.1:2*pi);
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            Z(i,j)=p.pendEnergy([X(i,j),Y(i,j)],p);
        end
    end
    
    for i = 1:length(t)
        e(i) = p.pendEnergy(q(i,:),p);
    end
    
    figure('units','normalized','outerposition',[0 0 1 1])
    for k = 1:length(t)
        hold on;
        surf(X,Y,Z,'FaceAlpha',0.2,'EdgeAlpha',0.1);
        plot3(q(1:k,1),q(1:k,2),e(1:k),'b','LineWidth',2)
        scatter3(q(k,1),q(k,2),e(k),100,'r','filled')
        plot(q(1:k,1),q(1:k,2),'b','LineWidth',2)
        scatter(q(k,1),q(k,2),100,'r','filled')
        contour(X,Y,Z,[1:25])
        view(20,20)
        drawnow;
        M(k) = getframe(gcf);
        clf;
    end
    %%
    
    video = VideoWriter('tmp.avi','Motion JPEG AVI');
    video.FrameRate = 24;
    open(video);
    writeVideo(video,M);
    close(video);
end