function [ym] = avgplot(x,y,f,f0,xlab,ylab,tit,power)

ym=max(y(:))/2;

clf
pcolor(x,y-ym,f)
set(gca,'Position',[0.1300    0.1100    0.550    0.8150])
title(tit,'fontsize',[14])
shading interp
colorbar('Location','WestOutside')
xlabel(xlab,'fontsize',[14])
ylabel(ylab,'fontsize',[14])
set(gca,'fontsize',[14])

axes('Position',[ 0.70    0.1100    0.23    0.815])
yavg=mean(y-ym,2);
favg=mean(f.^power,2).^(1./power);
f0avg=mean(f0.^power,2).^(1./power);

plot(favg,yavg,f0avg,yavg,'r--','LineWidth',2)
set(gca,'Position',[ 0.70    0.1100    0.23    0.815])
axis tight
set(gca,'YTickLabel',[])
set(gca,'fontsize',[14])

%title(tit,'fontsize',[14])

end