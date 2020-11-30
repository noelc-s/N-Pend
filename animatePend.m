%%
xline(1)
hold on;
h=plot(0,0,'MarkerSize',10,'Marker','.','LineWidth',2);
range=1.1*sum(p.l); axis([-range range -range range]); axis square;
if p.drawTail
    tail = [];
    tailPlot = plot(0,0,'LineWidth',2);
end
set(gca,'nextplot','replacechildren');

for i=1:length(q)-1
    Xcoord = zeros(p.n+1,1);
    Ycoord = zeros(p.n+1,1);
    for j = 1:p.n
        for k = 1:j
            Xcoord(j+1) = Xcoord(j+1)+p.l(k)*sin(q(i,k));
            Ycoord(j+1) = Ycoord(j+1)-p.l(k)*cos(q(i,k));
        end
    end
    set(h,'XData',Xcoord,'YData',Ycoord);
    tail = [tail; Xcoord(end) Ycoord(end)];
    set(tailPlot,'XData',tail(:,1),'YData',tail(:,2));
    axis off
    drawnow;
%     M(i) = getframe(gcf);
    pause(t(i+1)-t(i));
end
%%
% video = VideoWriter('tmp.avi','Uncompressed AVI');
% video.FrameRate = 24;
% open(video);
% writeVideo(video,M);
% close(video);