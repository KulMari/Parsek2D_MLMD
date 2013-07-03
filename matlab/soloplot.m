function [ym] = soloplot(x,y,f,xlab,ylab,tit)

ym=max(y(:))/2;

hold off
fmin1=min(f(:));
fmax1=max(f(:));
pcolor(x,y-ym,f)
title(tit,'fontsize',[14])
shading interp
colorbar
xlabel(xlab,'fontsize',[14])
ylabel(ylab,'fontsize',[14])

end