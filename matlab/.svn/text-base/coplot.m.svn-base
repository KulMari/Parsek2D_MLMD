function [ym] = coplot(x,y,f,ay,xlab,ylab,tit)

ym=max(y(:))/2;

hold off
fmin1=min(f(:));
fmax1=max(f(:));
pcolor(x,y-ym,f)
title(tit,'fontsize',[14])
shading interp
colorbar
hold on
aymin=min((ay(:)));
aymax=max((ay(:)));
fmax=max(fmax1,fmin1)
fmin=min(fmin1,fmax1)
ay2=((ay-aymin)/(aymax-aymin)*(fmax-fmin)+fmin);

contour(x,y-ym,ay2,30,'w')
xlabel(xlab,'fontsize',[14])
ylabel(ylab,'fontsize',[14])

end