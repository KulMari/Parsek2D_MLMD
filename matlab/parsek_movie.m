scrsz = get(0,'ScreenSize');
%figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])

numframes = 40;
num_frames_per_second = 2;
[nx_grid,ny_grid,time_slices] = size(Bx)
%A=moviein(time_slices); 
aviobj = avifile ( 'B_magnitudo_mms7.avi', 'fps', num_frames_per_second,'quality',100); 
set(gca,'NextPlot','replacechildren')

x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);
[X,Y] = meshgrid(x,y);
for i=1:time_slices
%cav = curl(X,Y,Bx(:,:,i),By(:,:,i));
%cav = Ex(:,:,i).^2 + Ey(:,:,i).^2;
%cav = rho(:,:,i);
cav = Bx(:,:,i).^2 + By(:,:,i).^2 + Bz(:,:,i).^2;
%contourf(X,Y,cav); 

shading interp
surf(X,Y,cav,'edgecolor','none','facecolor','blue')
caxis([4E-4 8.5E-4]);
%colormap('Jet');
axis([0 Lx 0 Ly 0 8E-3]);
lighting phong
shading interp
camlight(0,90)
%camlight(0,-90)
view(2)
%hold on;
%quiver(X,Y,Bx(:,:,i),By(:,:,i),'k-')     
%xlabel('x');
%ylabel('y');
title(i);
colormap('Jet');
colorbar
frame=getframe(gca);
aviobj = addframe ( aviobj, frame );
%hold off
%cav = curl(X,Y,Bx(:,:,i),By(:,:,i));
end 
aviobj = close ( aviobj );
%save BxBy.mat A
%figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
%movie2avi(A,'rek.avi');

% surf(cav,'edgecolor','none','facecolor','blue')
% 
% colormap('Hot');
% lighting phong
% shading interp
% camlight(0,90)
% camlight(0,-90)
% view(2)
%title('J intensity')
%xlabel('x')
%ylabel('y')

