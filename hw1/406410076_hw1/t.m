figure, subplot(2,2,1),plot(rand(3)),set(gca,'xtick',[],'ytick',[])
 subplot(2,2,2),plot(rand(3)),set(gca,'xtick',[],'ytick',[])
 subplot(2,2,3),plot(rand(3)),set(gca,'xtick',[],'ytick',[])
 subplot(2,2,4),plot(rand(3)),set(gca,'xtick',[],'ytick',[])
 ha=get(gcf,'children');
 set(ha(1),'position',[.5 .1 .4 .4])
 set(ha(2),'position',[.1 .1 .4 .4])
 set(ha(3),'position',[.5 .5 .4 .4])
 set(ha(4),'position',[.1 .5 .4 .4])