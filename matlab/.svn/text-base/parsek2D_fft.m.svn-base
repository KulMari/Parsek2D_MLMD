%% Fourier decomposition and growth rate computation

% quasi-1D geometry

Nmodes = 3; % number of modes to show (N = 0,...,Nmodes)
%EndDt = 60; % time in which modes saturate
EndCycle=150;
N=Nmodes+1;

indexes=[];
for i=0:Nmodes
    indexes=[indexes {num2str(i)}];
end

Bx0=mean(mean(Bx(:,:,1)));
By0=mean(mean(By(:,:,1)));

ref=sqrt(Bx0^2+By0^2);


Ex_fft=fft(squeeze(mean(Ex./ref)))';
Ey_fft=fft(squeeze(mean(Ey./ref)))';
Ez_fft=fft(squeeze(mean(Ez./ref)))';

By_fft=fft(squeeze(mean((By-By0)./ref)))';
Bz_fft=fft(squeeze(mean(Bz./ref)))';

E_fft=abs(Ez_fft);%+abs(Ey_fft)+abs(Ez_fft);
B_fft=abs(Bz_fft);%+abs(Bz_fft);

% Compute growth-rate for each mode
%EndCycle=ceil(EndDt/Dt);
clear growth
for i=1:N
    for k=10:EndCycle
        growth(i,k-9)=log(By_fft(k,i)/By_fft(1,i))/(k-1)/Dt/single(Bx_time(2))/0.098521; % note this is normalized with the same normalization of Dt
    end
end

% Compute omega(k)
Ex_fft2=fft2(squeeze(mean(Ex./ref))');
Ey_fft2=fft2(squeeze(mean(Ey./ref))');
Ez_fft2=fft2(squeeze(mean(Ez./ref))');
Ex_fft2=fftshift(Ex_fft2);
Ey_fft2=fftshift(Ey_fft2);
Ez_fft2=fftshift(Ez_fft2);

FirstM=size(E_fft,2)/2+1-Nmodes;
LastM=size(E_fft,2)/2+2+Nmodes;
FirstT=size(E_fft,1)/2+1-10;
LastT=size(E_fft,1)/2+2+10;

Times=single(Ex_time)*Dt*0.098521;
figure;
subplot(2,2,1)
semilogy(Times,B_fft(:,1:N));
legend(indexes,'Location','SouthEast');
xlabel('T');ylabel('\delta E / B_0')
title('Mode decomposition');

subplot(2,2,2)
pcolor(1:size(E_fft,2),Times,fftshift(E_fft,2));shading flat;colorbar
axis([ FirstM LastM 0 Times(end)])
set(gca,'XTick',[FirstM-1+0.5:LastM-0.5])
set(gca,'XTickLabel',[-N:N])
xlabel('Mode');ylabel('T')

subplot(2,2,3)
plot(Times(10:EndCycle),growth')
legend(indexes,'Location','NorthEast');
title('Growth rate')
xlabel('T');ylabel('\gamma')

subplot(2,2,4)
pcolor([FirstM:LastM],[FirstT:LastT],abs(Ey_fft2(FirstT:LastT,FirstM:LastM)+(Ez_fft2(FirstT:LastT,FirstM:LastM))+(Ex_fft2(FirstT:LastT,FirstM:LastM))));colorbar
set(gca,'XTick',[FirstM-1+0.5:LastM-0.5])
set(gca,'XTickLabel',[-N:N])
set(gca,'YTick',[FirstT-1+0.5:2:LastT-0.5])
set(gca,'YTickLabel',[-11:2:11])
title('Dispersion relation \omega(k) in Fourier space');
xlabel('k');ylabel('\omega')
