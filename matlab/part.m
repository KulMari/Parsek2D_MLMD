addpath 'matlab-enrico'


close all

global results_dir variable_list

if(0==1)
    clear all
    results_dir='/home/gianni/cluster/storage/parsek/mms/results/';
    %results_dir='/home/gianni/franklin/results/';
    %results_dir='/Users/gianni/Documents/runni/parsek/results/';
    variable_list='x v ID B J rho';

    parsek2D
end

elettroni=0

[np ntp] = size(x0);
[nx ny nt]=size(Bx);

wci=.0195;dt=.25;
vthe=.045
va=.0195
time=x_time*wci*dt;
aver=0;

for itp=1:ntp
    
        it=1+(itp-1)*nt/ntp

    
    bbx=squeeze(sum(Bx(1:end,1:end,it-aver:it),3))/(aver+1);
bby=squeeze(sum(By(1:end,1:end,it-aver:it),3))/(aver+1);
bbz=squeeze(sum(Bz(1:end,1:end,it-aver:it),3))/(aver+1);

if(elettroni)
jsx0=squeeze(sum(Jxs0(1:end,1:end,it-aver:it),3))/(aver+1);
jsy0=squeeze(sum(Jys0(1:end,1:end,it-aver:it),3))/(aver+1);
jsz0=squeeze(sum(Jzs0(1:end,1:end,it-aver:it),3))/(aver+1);
rrho=squeeze(sum(rhos0(1:end,1:end,it-aver:it),3))/(aver+1);
else
jsx0=squeeze(sum(Jxs1(1:end,1:end,it-aver:it),3))/(aver+1);
jsy0=squeeze(sum(Jys1(1:end,1:end,it-aver:it),3))/(aver+1);
jsz0=squeeze(sum(Jzs1(1:end,1:end,it-aver:it),3))/(aver+1);
rrho=squeeze(sum(rhos1(1:end,1:end,it-aver:it),3))/(aver+1);   
end
bbb=sqrt(bbx.^2+bby.^2+bbz.^2);
upar=(jsx0.*bbx+jsy0.*bby+jsz0.*bbz)./(rrho.*bbb);

[xx yy]=meshgrid(1:ny,1:nx);

ay=vecpot(xx,yy,bbx,bby);
[aymin,imin]=min(ay(nx/2,:))
xx=xx/ny*20;
yy=yy/nx*10;

caso=2;
if(caso==1)
xzero=.5*(xx(nx/2,imin)+xx(nx/2,imin+1))
yzero=5;
else
    xzero=13.25;
    yzero=5-1.4;
end
    xmin=xzero-.5;
    xmax=xzero+.5;
    ymin=yzero-.5;
    ymax=yzero+.5;
    
h=figure(3)
set(h,'Position' , [565 50 560 420]);
coplot(xx,yy,upar/va,ay,'x/d_i','z/d_i',['U_{e||}/v_a(\omega_{ci}t=' num2str(time(itp)) ')'])
ym=max(yy(:))/2;
plot([xmin xmax xmax xmin xmin],[ymin ymin ymax ymax ymin]-ym,'m','linewidth',[3])

print('-dpng','-r300',['fgE' num2str(itp) '.png'])
if (elettroni) 
    ii=x0(:,itp)>xmin &x0(:,itp)<xmax &y0(:,itp)<ymax &y0(:,itp)>ymin ;
    wpsub=w0(ii,itp);
else
    ii=x1(:,itp)>xmin &x1(:,itp)<xmax &y1(:,itp)<ymax &y1(:,itp)>ymin ;
    wpsub=w1(ii,itp);
end    
    h=figure(1)
    set(h,'Position' , [5 50 560 420]);
    hist(wpsub/va,150)
    xlabel('v_p/v_a','fontsize',[14])
    ylabel('counts','fontsize',[14])
  
    title(['\omega_{ci}t=' num2str(time(itp))],'fontsize',[14])
    
    print('-dpng','-r300',['fgF' num2str(itp) '.png'])

pause(.1)
end    