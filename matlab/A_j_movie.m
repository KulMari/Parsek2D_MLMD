scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])


[nx_grid,ny_grid,time_slices] = size(Bx)
A=moviein(time_slices); 
set(gca,'NextPlot','replacechildren')

x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);
[X,Y] = meshgrid(x,y);
for i=1:time_slices
%cav = curl(X,Y,Bx(:,:,i),By(:,:,i));
cav = Bz(:,:,i);
contourf(X,Y,cav); 

%caxis([-0.07 0]);
%shading interp
hold on;
quiver(X,Y,Bx(:,:,i),By(:,:,i))     
xlabel('x');
ylabel('y');
title(i);
colormap('Hot');
colorbar
A(:,i)=getframe;
hold off
end 

figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
%movie2avi(A,'rek.avi');
surf(-cav,'edgecolor','none','facecolor','blue')

colormap('Hot');
lighting phong
shading interp
camlight(0,90)
camlight(0,-90)
view(2)
title('J intensity')
xlabel('x')
ylabel('y')