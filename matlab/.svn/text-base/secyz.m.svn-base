addpath 'matlab-enrico'


close all

global results_dir variable_list

if(0==0)
    clear all
    %results_dir='/home/gianni/cluster/storage/parsek/mms/results/';
    %results_dir='/home/gianni/cluster/storage/parsek/mms/prit/';
    %results_dir='/home/gianni/cluster/storage/parsek/selffeding/results/';
    %results_dir='/home/gianni/cluster/storage/parsek/selffeding/prit/';

    %results_dir='/home/gianni/franklin/results/';
    results_dir='/home/gianni/Documents/workspace/Parsek2D/results/';
    variable_list='E B rho J';

    parsek2D
end

wci=1
L=1

[nx ny nt]=size(Bx);

it=nt

iy=2;

% for mms
%wci=.0195;

%for self-feeding
%wci=.069;

wci=max([Bx0 By0 Bz0])
% for LH frequency
wci=wci*sqrt(abs(qom(1)));
% for electron plasma frequency 
%wci=1;

vthe=uth(2);
va=wci;

time=double(Bx_time)*wci*Dt;
aver=0;

recon=[];

indexf=1
first=1;
for it=1:nt

   
bbx=squeeze(sum(Bx(1:end,1:end,it-aver:it),3))/(aver+1);
bby=squeeze(sum(By(1:end,1:end,it-aver:it),3))/(aver+1);
bbz=squeeze(sum(Bz(1:end,1:end,it-aver:it),3))/(aver+1);
eex=squeeze(sum(Ex(1:end,1:end,it-aver:it),3))/(aver+1);
eey=squeeze(sum(Ey(1:end,1:end,it-aver:it),3))/(aver+1);
eez=squeeze(sum(Ez(1:end,1:end,it-aver:it),3))/(aver+1);
jsx0=squeeze(sum(Jxs0(1:end,1:end,it-aver:it),3))/(aver+1);
jsx1=squeeze(sum(Jxs0(1:end,1:end,it-aver:it),3))/(aver+1);
jsy0=squeeze(sum(Jys0(1:end,1:end,it-aver:it),3))/(aver+1);
jsz0=squeeze(sum(Jzs0(1:end,1:end,it-aver:it),3))/(aver+1);
rrho=squeeze(sum(rhos0(1:end,1:end,it-aver:it),3))/(aver+1);
rrhoi=squeeze(sum(rhos1(1:end,1:end,it-aver:it),3))/(aver+1);

bbb=sqrt(bbx.^2+bby.^2+bbz.^2);
upar=(jsx0.*bbx+jsy0.*bby+jsz0.*bbz)./(rrho.*bbb);
epar=(eex.*bbx+eey.*bby+eez.*bbz)./(bbb);

[xx yy]=meshgrid(1:ny,1:nx);

ay=vecpot(xx,yy,bbx,bby);
xx=xx/ny*Lx;
yy=yy/nx*Ly;

h=figure(1);
set(h,'Position' , [5 5 560 420]);
var=jsx0./rrho;
if(first==1) 
    var0=var;
end
avgplot(xx,yy,var,var0,'x/d_i','y/d_i',['U_{xe}(\omega_{lh}t=' num2str(time(it),'%10.2f') ')'],1)
F1(indexf)=getframe(gcf);

h=figure(2)
set(h,'Position' , [5 550 560 420]);
var=rrho*4*pi;
if(first==1) 
    var1=var;
end
avgplot(xx,yy,var,var1,'x/d_i','y/d_i',['\rho_{e}(\omega_{lh}t=' num2str(time(it),'%10.2f') ')'],1)
F2(indexf)=getframe(gcf);

h=figure(3)
set(h,'Position' , [565 5 560 420]);
var=bbz;
if(first==1) 
    var2=var;
end
var=jsx1./rrhoi;
if(first==1) 
    var2=var;
end
avgplot(xx,yy,var,var2,'x/d_i','y/d_i',['U_{xi}(\omega_{lh}t=' num2str(time(it),'%10.2f') ')'],1)
F3(indexf)=getframe(gcf);

h=figure(4)
set(h,'Position' , [565 550 560 420]);
var=eey;
if(first==1) 
    var3=var;
end
avgplot(xx,yy,var,var3,'x/d_i','y/d_i',['E_{y}(\omega_{lh}t=' num2str(time(it),'%10.2f') ')'],2)
F4(indexf)=getframe(gcf);

h=figure(5)
set(h,'Position' , [1125 550 560 420]);
var=rrhoi*4*pi;
if(first==1) 
    var4=var;
end
avgplot(xx,yy,var,var4,'x/d_i','y/d_i',['\rho_{i}(\omega_{lh}t=' num2str(time(it),'%10.2f') ')'],1)
F5(indexf)=getframe(gcf);

% surf(xx,yy-max(yy(:))/2,eez,'edgecolor','none','facecolor','blue')
% lighting phong
% shading interp
% camlight(0,90) % luce dall'alto
% view(2) %visione dall'alto
% axis tight

indexf=indexf+1;

recon=[recon;max(ay(end/2,:))-min(ay(end/2,:))];

h=figure(100)
set(h,'Position' , [1125 5 560 420]);
plot(recon/va)
xlabel('\omega_{ci}t','fontsize',[14])
ylabel('\Delta \Psi/B_0d_i','fontsize',[14])
pause(.1)

first=0;

end

return

name='film_mms_back'
movie2avi(F1,[name num2str(1) '.avi'],'fps',[2],'quality',[100])
%!ffmpeg -i film_mms1.avi -sameq -r 24 film_mms1.mpg
%!ffplay film_mms1.mpg
movie2avi(F2,[name num2str(2) '.avi'],'fps',[2],'quality',[100])
movie2avi(F3,[name num2str(3) '.avi'],'fps',[2],'quality',[100])
movie2avi(F4,[name num2str(4) '.avi'],'fps',[2],'quality',[100])
movie2avi(F5,[name num2str(5) '.avi'],'fps',[2],'quality',[100])
!ffmpeg -i film_mms_back1.avi -sameq -r 24 film_mms_backe1.mpg
!ffmpeg -i film_mms_back2.avi -sameq -r 24 film_mms_backe2.mpg
!ffmpeg -i film_mms_back3.avi -sameq -r 24 film_mms_backe3.mpg
!ffmpeg -i film_mms_back4.avi -sameq -r 24 film_mms_backe4.mpg
!ffmpeg -i film_mms_back5.avi -sameq -r 24 film_mms_backe5.mpg