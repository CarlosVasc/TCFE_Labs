close all
clear all
clc


%------------------------
f = logspace(0,8,800);

C1 = 220e-9
C2 = 110e-9
R1 = 1e3
R2 = 1e3
R3 = 1e3
R4 = 100e3


gain = (R1.*2.*pi.*f.*C1./(1+R1.*j.*2.*pi.*f.*C1)).*(1 + R4/R3).*(1./(j.*2.*pi.*f.*C2.*R2 + 1));
gain_db=20*log10(abs(gain));

I_max=find (gain_db==max(gain_db));
f_max=f(I_max);

gain_max=max(gain);
gainDB_max=max(gain_db);
Z_in = (R1 + 1/(j*2*pi*f_max*C1))
Z_out = R2/(j*2*pi*f_max*C2)/(R2+1/(j*2*pi*f_max*C2))

Zin_a=abs(Z_in)
Zout_a=abs(Z_out)

%---------------------
G_desvio =abs(100- max(abs(gain)));
F_desvio=abs(1000-f_max);
cost_amp = ((8.661e-12 + 30e-12)*1e6 + (100000 + 5305 + 5305 + 1836 + 1836 + 13190000 + 50 + 100 + 18160)/1000 + (2*0.1)) + (0.220*2)
cost=cost_amp+(1e-3*(R1 + R2 + R3 + R4)+1e6*(C1 + C2))
Merit = 1/(cost * (G_desvio + F_desvio+1e-6))
%------------- Export ------------
fid = figure ();
semilogx(f, gain_db);
title("Gain");
legend("v_o(f)/v_i(f)");
ylabel("Gian (dB)")
xlabel("frequency (Hz)")
print(fid ,"Gain_teo.eps","-depsc");


fid2=fopen("valores.tex","w+");
fprintf(fid2,"C1 & %e\\\\ \\hline \n\
C2 & %.5e\\\\ \\hline \n\
R1 & %.5e\\\\ \\hline \n\
R2 & %.5e\\\\ \\hline \n\
R3 & %.5e\\\\ \\hline \n\
R4 & %.5e\\\\ \\hline", C1,C2,R1,R2,R3,R4)
fclose(fid2);


fid3=fopen("merit.tex","w+");
fprintf(fid3,"Gain Deviation & %f\\\\ \\hline \n\
Frequency deviation & %.f\\\\ \\hline \n\
Cost & %.5e\\\\ \\hline \n\
Merit & %.5e\\\\ \\hline", G_desvio,F_desvio,cost,Merit)
fclose(fid3);


fid4=fopen("Val_meio.tex","w+");
fprintf(fid4,"Gain & %.5e\\\\ \\hline \n\
Gain (db) & %.5e\\\\ \\hline \n\
Input impedance &  %.5e + %.5ej\\\\ \\hline \n\
Output impedance &  %e + %ej\\\\ \\hline", gain_max,gainDB_max,real(Z_in), imag(Z_in),real(Z_out) , imag(Z_out))
fclose(fid4);


