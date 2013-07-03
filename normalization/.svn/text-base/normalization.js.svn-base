function calculate(form) {
  var n0 =   form.text1.value;
  var B0 =  form.text2.value;
  var T  =  form.text3.value;
  var miome= form.text4.value;
  
  
  // formula
  
  var me=9.1094e-28;
  var mi=me*miome;  
  var kb=1.3807e-16;
  var q=4.8032e-10;
  var c=3e10;

  //  Quantities calculated

  T=T/1.1604e4; // T in eV;
  
  var Va=B0/Math.sqrt(4*Math.PI*n0*mi); // Alfven velocity
  
  var ld=7.43e2*Math.sqrt(T/n0);  // Debye length - cm
  var om_pi=Math.sqrt(4*Math.PI*n0*Math.pow(q,2)/mi); // ion plasma frequency - rad/s
  var om_pe=om_pi*Math.sqrt(miome); // electron plasma frequency - rad/s
  var om_ci=q*B0/mi/c; // ion cyclotron frequency - rad/s
  var om_ce=om_ci*miome; // electron cyclotron frequency - rad/s
  var vthe=4.19e7*Math.sqrt(2*T); 
  var vthi=vthe/Math.sqrt(miome);
  
  // wave stuff
  var Nk=form.text5.value/Math.cos(form.text6.value*(Math.PI/180));
  var k=Nk*om_pi/c;   // reference wavevector k
  var Nw=1;           // The box will be Nw wavelength long
  var L=2*Math.PI*Nw/k;    // Length of the box         

  var Ncell=form.text7.value;
 
  // Time and number of cycles;
  // The total time of simulation is given by N/growth-rate of a given
  // perturbation

  // Assign here N and the growth-rate gamma 
  var N=form.text8.value;
  var gamma=form.text9.value*om_ci;
  var Time=2*Math.PI*N/gamma;
  
  // Assign here how many time-steps for an electron gyromotion;
  var steps=form.text10.value ;
  var dt=1/steps/om_ce;

  var Ncycles=Time/dt;
  

  // Assign :
  var Npcel= Math.pow(10,2); // Number of particles per cell;
  var Nparticles=Math.pow(Npcel*Ncell,2); //
  var dx=L/Ncell;
  // Decide normalization quantities
  var Norm_time = "1/&omega<sub>pi</sub>";  
  var Norm_velocity = "c";
  var Norm_space = "&omega<sub>pi</sub>";


  var norm_time=1/om_pi;
  var norm_velocity= c;
  var norm_space= c/om_pi;
  var norm_qom=q/mi;  
  var norm_B0=1/(c/om_pi)*Math.pow(c,2)/(q/mi); 
  var norm_rho=norm_B0/(c/om_pi);
  

  var DT=dt/norm_time;
  var DX=dx/norm_space;
  var L_normalized=L/norm_space;
  var Time_normalized=Time/norm_time;
  var qom_i_normalized=q/mi/norm_qom;
  var qom_e_normalized=-q/me/norm_qom;
  var B0_normalized=B0/norm_B0;
  var rho=q*n0;
  var rho_normalized=rho/norm_rho;
  var vthi_normalized=vthi/norm_velocity;
  var vthe_normalized=vthe/norm_velocity;
  var c_normalized=c/norm_velocity;
  var Va_normalized=Va/norm_velocity;
  var om_pi_normalized=om_pi*norm_time;
  var om_pe_normalized=om_pe*norm_time;
  var om_ci_normalized=om_ci*norm_time;
  var om_ce_normalized=om_ce*norm_time;

  var k_normalized=k*norm_space;
  var delta_E=Math.pow(om_ce_normalized,1.5)/Math.sqrt(vthe_normalized);
  var delta_rho= Math.pow(om_ce_normalized,2.5)/Math.pow(vthe_normalized,1.5);
  var delta_j=Math.pow(om_ce_normalized,2.5)/Math.sqrt(vthe_normalized);

  // write to display
  document.write("<h2> Results </h2>");
  document.write("<hr>")
  document.write("<h3> Length scale</h3>");
  document.write("Physical size of the box is " + (L/1e5) + " km <br>");
  document.write("Debye length is  " + (ld/1e2) + " m <br>");
  document.write("Number of cell (dx = " +(dx/ld) +" l<sub>d</sub>) is " + Ncell + "<br>");
  document.write("<hr>")
  document.write("<h3> Time Scale</h3>");
  document.write("Total physical time of the simulation is " + Time +  " sec <br>");
  document.write("Physical dt is " + dt + " sec <br>");
  document.write("Number of cycles is " + (Math.ceil(Ncycles)) + "<br>");
  document.write("<hr>");
  document.write("<h3> Check simulation constraints</h3>");
  
  if ((vthe*dt/dx) > 1){
     document.write("v<sub>the</sub> &lt dx\dt  is <b>not satisfied</b> --> increase steps or Nw , or decrease Ncell or Nk <br>");
    } else {
    document.write("v<sub>the</sub> &lt dx\dt  is <b> satisfied</b> <br>");
  }

  if ((vthe*dt/dx)<0.1) {
    document.write("v<sub>the</sub> *dt/dx &gt 0.1  is <b>not satisfied</b> --> finite grid instability<br>");
 } else {
   document.write("v<sub>the</sub> *dt/dx &gt 0.1  is <b>satisfied</b><br>");
  }
  document.write("<hr>");
  // NORMALIZATION VALUES FOR 
  document.write("<h3> Normalization parameters</h3>");
  document.write("Velocity--> " + Norm_velocity + "<br>");
  document.write("Time --> " + Norm_time+ "<br>");
  document.write("qom<sub>i</sub> = " + qom_i_normalized + "<br>");
  document.write("qom<sub>e</sub> = " + qom_e_normalized + "<br>");
  document.write("L = " + L_normalized + "<br>");
  document.write("dx = " + DX + "<br>");
  document.write("n<sub>cell</sub> = " + Ncell + "<br>");
  document.write("Time = " + Time_normalized + "<br>");
  document.write("dt = " + DT+ "<br>");
  document.write("n<sub>cycles</sub> = " + Ncycles+ "<br>");
  document.write("B<sub>0</sub> = " + B0_normalized+ "<br>");
  document.write("&rho = " + rho_normalized + "<br>");
  document.write("c = " + c_normalized + "<br>");
  document.write("V<sub>a</sub> = " + Va_normalized + "<br>");
  document.write("v<sub>thi</sub> = " + vthi_normalized + "<br>");
  document.write("v<sub>the</sub> = " + vthe_normalized + "<br>");
  document.write("&omega<sub>pi</sub> = " + om_pi_normalized + "<br>");
  document.write("&omega<sub>pe</sub> = " + om_pe_normalized + "<br>");
  document.write("&omega<sub>ci</sub> = " + om_ci_normalized + "<br>");
  document.write("&omega<sub>ce</sub> = " + om_ce_normalized + "<br>");
  document.write("&Delta E = " + delta_E + "<br>");
  document.write("&Delta &rho = " + delta_rho + "<br>");
  document.write("&Delta j = " + delta_j + "<br>");
  document.write("1 electron gyroradius is " + (0.5*vthe_normalized/om_ce_normalized/DX) + " cell size <br>");
  document.write("v<sub>the</sub>*dt/dx = " + (vthe*dt/dx) + "<br>");
  document.write("<hr>");
  
}

