close all
clear all
pkg load symbolic
pkg load control
pkg load signal

Vfonte = 230;
f=50;
w=2*pi*f;

R1 = 20000;
C = 10e-6;
R2 = 2000;
voltas = 17;

A = Vfonte/voltas;
t=linspace(0, 0.2, 1000);

vs=A*cos(w*t); %depois de transformador
toff = 1/w * atan(1/(w*R1*C));
v_fex = A*cos(w*toff)*exp(-(t-toff)/(R1*C));

for i=1:length(t)
	  v0hr(i) = abs(vs(i));
end

for i=1:length(t)
  if t(i) < toff
    v0(i) = v0hr(i);
  elseif v_fex(i) > v0hr(i)
    v0(i) = v_fex(i);
  else 
    toff = toff + 1/(2*f) ;
    v_fex = A*abs(cos(w*toff))*exp(-(t-toff)/(R1*C));
    v0(i) = v0hr(i);
  end
end

average = mean(v0);
ripple = max(v0) - min(v0);


%Regulador

n_diodes = 17;
Von = 0.7;
v0_2_dc = Von*n_diodes;
vt = 0.025;
Is = 1e-14;
n = 1;
Rd = n*vt/(Is*exp(Von/(n*vt)));

for i = 1:length(t)
    v0_2_ac(i) = n_diodes*Rd/(n_diodes*Rd+R2) * (v0(i)-average);
end

v0_2 = v0_2_dc + v0_2_ac;
averageR = mean(v0_2);
rippleR = max(v0_2)-min(v0_2); 
cost = R1/1000 + R2/1000 + (n_diodes+0.4)*0.1+ C*1e6 ; 
merito = 1/(cost*(rippleR + abs(averageR - 12) + 1e-6))


%grafs
fid1 = figure();
plot (t*1000, vs, ";vs in transformer;", t*1000,v0, ";v in envelope;", t*1000,v0_2, ";v in regulator;");
xlabel ("t(ms)")
ylabel ("v (Volts)")
legend('Location','southeast');
print (fid1, "V_todas.eps", "-depsc");


fid2 = figure();
plot (t*1000,v0_2-12);
xlabel ("t(ms)")
ylabel ("v (Volts)")
print (fid2, "diferenca.eps", "-depsc");

