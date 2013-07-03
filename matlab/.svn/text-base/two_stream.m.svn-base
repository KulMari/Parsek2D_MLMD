scrsz = get(0,'ScreenSize');
%figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
num_frames_per_second = 8;
%A=moviein(time_slices); 
aviobj = avifile ( 'v05_phase_space.avi', 'fps', num_frames_per_second ); 
set(gca,'NextPlot','replacechildren')

Lx = 12.0
[np,time] = size(x0)
% momentum
p = zeros(np,1);
v = zeros(np,1);
% gamma 
gamma = zeros(np,1);
%define the BOX
vYEdge = linspace(0,Lx,128);
vXEdge = linspace(-1.8,1.8,128);

for i=1:time-1
  i
%   plot(-rhos0(:,i),'k','LineWidth',1.0);
%   hold on
%   plot(rhos1(:,i),'r','LineWidth',1.0);
%   axis([0 512 0 0.2]);
%   hold off
  % calculate the classical momentum
  p = (u0(:,i));
  % calculate the relativistic gamma. In the PARSEK normalization
  gamma = 1./sqrt(1-u0(:,i).^2);
  p = gamma.*p;
  max(p(:,1))
  min(p(:,1))
%   mHist2d = hist2d([x0(:,i), p(:,1)],vYEdge,vXEdge);
  
  vector1 = find(p>0);
  vector2 = find(p<0);
  
  plot(x0(vector1,i),p(vector1,1),'k.','MarkerSize',.6);
  hold on
  plot(x0(vector2,i),p(vector2,1),'k.','MarkerSize',.6);
  hold off;
 % hold on;
%   % calculate the classical momentum
   %p = (q1(:,i).*u1(:,i));
%   % calculate the relativistic gamma. In the PARSEK normalization
   %gamma = 1./sqrt(1-u1(:,i).^2);
%   max(p(:,1))
  
%    plot(x1(:,i),p(:,1),'r.','MarkerSize',.1);
%    hold off;
 whitebg([0.8 0.8 0.8]);
 axis([0 Lx -1.8 1.8]);
  %axis([0 Lx -2E-6 2E-6]);
%   surf(mHist2d','edgecolor','none','facecolor','blue')
%   axis([0 128 0 128]);
%    colorbar;
%    caxis([0 30]);
%   lighting phong
%   shading interp
  camlight(0,90)
  view(2)
  colormap('Jet');


frame=getframe(gca);
aviobj = addframe ( aviobj, frame );
end

aviobj = close ( aviobj );
