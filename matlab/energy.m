t=single(k_nrg_time)*Dt*0.098521;
Dvol=Dx*Dy;

total=(E_energy+B_energy)*Dvol+k_energy_total;

figure
subplot(2,1,1)
plot(t,k_energy0./k_energy0(1))
grid;hold on
plot(t,k_energy1./k_energy1(1),'k')
plot(t,k_energy_total./k_energy_total(1),'r');
plot(t,total./total(1),'m')
legend('electron','ion','total K','Total')
xlabel('T\Omega_p');
ylabel('arb')

subplot(2,1,2)
plot(t,B_energy*Dvol./total)
grid;hold on
%plot(t,(B_energy*Dvol+k_energy_total)./total,'r')

