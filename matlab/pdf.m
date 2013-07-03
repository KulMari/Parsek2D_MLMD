%% 2D contour plot of distribution function
X = Lx/2; Y = Ly/2; % point where to center the pdf
%X=6;Y=1.7;
box_Lx= Lx/400 ; box_Ly= Ly/3; % cannot be > Lx, Ly (dimension of box)
box_Lx= Lx ; box_Ly= Ly; 
s=0; % species in consideration
T=1; % time
sample =1 ; % consider one every 'sample' particles (ie sample =1 --> all particles)

%%%%%%%%%%%%%%%%%%%%%%%
u=eval(['u' num2str(s)]);
v=eval(['v' num2str(s)]);
w=eval(['w' num2str(s)]);
x=eval(['x' num2str(s)]);
y=eval(['y' num2str(s)]);
q=eval(['q' num2str(s)]);

list=[];

% collecting particles on which compute pdf
if (box_Lx ~= Lx) || (box_Ly ~= Ly)
 for i=1:sample:size(u,1)
  if (X-0.5*box_Lx <= x(i,T)) && (x(i,T) <= X+0.5*box_Lx) && (Y-0.5*box_Ly <= y(i,T)) &&( y(i,T) <= Y+0.5*box_Ly)
   list=[list; q(i,T) u(i,T) v(i,T) w(i,T) sqrt(v(i,T)^2+w(i,T)^2)]; 
   
  end
 end
else
%    list=[q(:,T) u(:,T) v(:,T) w(:,T) sqrt(v(:,T).^2+w(:,T).^2)];
    list=[q(:,T) u(:,T)*cosd(80)-v(:,T)*sind(80) u(:,T)*sind(80)+v(:,T)*cosd(80) w(:,T) sqrt((u(:,T)*sind(80)+v(:,T)*cosd(80)).^2+w(:,T).^2)];
%    list=[q(:,T) u(:,T) q(:,T) w(:,T) sqrt(v(:,T).^2+w(:,T).^2)];

end
Nparticles=size(list,1);
disp(['Distribution function computed on ' num2str(Nparticles) ' particles'])
maxu=max(list(:,2)); minu=min(list(:,2));maxu=max(maxu,-minu);minu=-maxu;
maxv=max(list(:,3)); minv=min(list(:,3));maxv=max(maxv,-minv);minv=-maxv;
maxw=max(list(:,4)); minw=min(list(:,4));maxw=max(maxw,-minw);minw=-maxw;
maxvperp=max(list(:,5)); minvperp=0;

U=linspace(minu,maxu,30);dU=U(2)-U(1);
V=linspace(minv,maxv,30);dV=V(2)-V(1);
W=linspace(minw,maxw,30);dW=W(2)-W(1);
Vperp=linspace(minvperp,maxvperp,30);dVperp=Vperp(2)-Vperp(1);

XY=zeros(size(U,2),size(V,2));
XZ=zeros(size(U,2),size(W,2));
ZY=zeros(size(W,2),size(V,2));
XYZ=zeros(size(U,2),size(Vperp,2));

for i=1:Nparticles
if (~isnan(list(i,1)))
XY(floor(abs(list(i,2)-minu-1e-10)/dU)+1,floor(abs(list(i,3)-minv-1e-10)/dV)+1)=XY(floor(abs(list(i,2)-minu-1e-10)/dU)+1,floor(abs(list(i,3)-minv-1e-10)/dV)+1)+abs(list(i,1));
XZ(floor(abs(list(i,2)-minu-1e-10)/dU)+1,floor(abs(list(i,4)-minw-1e-10)/dW)+1)=XZ(floor(abs(list(i,2)-minu-1e-10)/dU)+1,floor(abs(list(i,4)-minw-1e-10)/dW)+1)+abs(list(i,1));
ZY(floor(abs(list(i,4)-minw-1e-10)/dW)+1,floor(abs(list(i,3)-minv-1e-10)/dV)+1)=ZY(floor(abs(list(i,4)-minw-1e-10)/dW)+1,floor(abs(list(i,3)-minv-1e-10)/dV)+1)+abs(list(i,1));
XYZ(floor(abs(list(i,2)-minu-1e-10)/dU)+1,floor(abs(list(i,5)-minvperp-1e-10)/dVperp)+1)=XYZ(floor(abs(list(i,2)-minu-1e-10)/dU)+1,floor(abs(list(i,5)-minvperp-1e-10)/dVperp)+1)+abs(list(i,1));
end

end
figure
%figure(1)
subplot(2,2,1)
contour(U+dU/2,V+dV/2,XY'./max(max(XY)));colorbar;axis equal
title('pdf (v_x,v_y)');xlabel('v_x');ylabel('v_y');grid on

%figure(2)
subplot(2,2,2)
contour(U+dU/2,W+dW/2,XZ'./max(max(XZ)));colorbar;axis equal
title('pdf (v_x,v_z)');xlabel('v_x');ylabel('v_z');grid on

%figure(3)
subplot(2,2,3)
contour(W+dW/2,V+dV/2,ZY'./max(max(ZY)));colorbar;axis equal
title('pdf (v_z,v_y)');xlabel('v_z');ylabel('v_y');grid on

subplot(2,2,4)
contour(U+dU/2,Vperp+dV,XYZ'./max(max(XYZ)));colorbar;axis equal
title('pdf (v_z,v_y)');xlabel('v_{||}');ylabel('v_{\perp}');grid on
