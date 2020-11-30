%%
figure('units','normalized','outerposition',[0 0 .5 .8])
range=1.1*p.n*p.l; axis([-range range -range range]); axis square;
set(gca,'nextplot','replacechildren');
hold on;
for z = 1:2
    h{z}=plot(0,0,'MarkerSize',20,'Marker','.','LineWidth',3);
end
axis off

for i=1:length(state)-1
    for z = 1:2
    Xcoord = zeros(n+1,1);
    Ycoord = zeros(n+1,1);
    for j = 1:n
        for k = 1:j
            Xcoord(j+1) = Xcoord(j+1)+p.l*sin(state(i,z*15,k));
            Ycoord(j+1) = Ycoord(j+1)-p.l*cos(state(i,z*15,k));
        end
    end
    set(h{z},'XData',Xcoord,'YData',Ycoord);
    drawnow;
    M(i) = getframe(gcf);
%     pause(t(i+1)-t(i));
end
end
%%
% video = VideoWriter('ic.avi','Uncompressed AVI');
% video.FrameRate = 60;
% open(video);
% writeVideo(video,M);
% close(video);