close all
clear all

pkg load symbolic
pkg load control
pkg load signal

%tirar dados do ficheiro.txt produzido pelo make
format longg
%abertura do ficheiro
fid=fopen("testi.txt","r");
ll=1;
while ~feof(fid)
     linha=fgetl(fid);
   if(linha)~=0 
       if ll>6 
        [dado,resto]=strtok(linha," ");
        [lixo,valor]=strtok(resto," ");
        d(ll)=str2double(valor);
        
       elseif ll==6
        [dado,resto]=strtok(linha," R1 = ");
        [lixoaux,valoraux]=strtok(resto," ");  
        [lixo,valor]=strtok(valoraux," ");
        d(ll)=str2double(valor);
        
       end
   else 
       d(ll)=0; 
   end
   
    ll=ll+1;
end
fclose(fid);

R1=d(6)*1000;
R2=d(7)*1000;
R3=d(8)*1000;
R4=d(9)*1000;
R5=d(10)*1000;
R6=d(11)*1000;
R7=d(12)*1000;
Vs=d(13);
C=d(14)*0.000001;
Kb=d(15)*0.001;
Kd=d(16)*1000;

G1=1/R1;
G2=1/R2;
G3=1/R3;
G4=1/R4;
G5=1/R5;
G6=1/R6;
G7=1/R7;

syms t

%----------  1  ----------------
A = [    0   0   1   0   0   0   0   0   0   0   0   0;
        -1   0   0   1   0   -1  0   0   0   0   0   0; 
        -Kb   0   0   0   0   0   0   0   0   1  0   0;
        0   -1   0   0   0   0   0   0   0   0   0   Kd;
        0   -1   0   0   0   1   0   0  -1   0   0   0;
        0    0   0   0   0   0   0 -G6   0   0   0  -1;
        0    0   0   0   0   0   0   0   0   0   1   0;
        0    0  -G1 G1+G2+G3   -G2   -G3   0   0       0   0    0   0;
        0   0   0  -G2-Kb   G2   Kb      0     0       0   0   0   0;
        0   0   0    Kb     0   -Kb-G5   G5    0       0   0   0   0;
        0   0   0    0      0    0       0   G7+G6   -G7   0   0   0; 
        0   0   G1  -G1     0   -G4      0     0       0   0   0   1];
    
B = [Vs;0;0;0;0;0;0;0;0;0;0;0];

NC=linsolve(A,B);
 

Vb= NC(1);
Vd= NC(2);
V1= NC(3);
V2= NC(4);
V3= NC(5);
V5= NC(6);
V6= NC(7);
V7= NC(8);
V8= NC(9);
Ib= NC(10);
Ic= NC(11);
Id= NC(12);

I_R1 = G1*(V1 - V2)
I_R2 = G2*(V2 - V3)
I_R3 = G3*(V5 - V2)
I_R4 = V5*G4
I_R5 = G5*(V6 - V5)
I_R6 = -V7*G6
I_R7 = G7*(V7 - V8)


%Criar ficheiro octave teorico
fid11=fopen("teorico1.tex","w+");
fprintf(fid11,"V1 & %f\\\\ \\hline \n\
V2 & %f\\\\ \\hline \n\
V3 & %f\\\\ \\hline \n\
V5 & %f\\\\ \\hline \n\
V6 & %f\\\\ \\hline \n\
V7 & %f\\\\ \\hline \n\
V8 & %f\\\\ \\hline \n\
IR1 & %f\\\\ \\hline \n\
IR2 & %f\\\\ \\hline \n\
IR3 & %f\\\\ \\hline \n\
IR4 & %f\\\\ \\hline \n\
IR5 & %f\\\\ \\hline \n\
IR6 & %f\\\\ \\hline \n\
IR7 & %f\\\\ \\hline", V1,V2,V3,V5,V6,V7,V8,I_R1,I_R2,I_R3,I_R4,I_R5,I_R6,I_R7)
fclose(fid11)


fid12=fopen("ngspice_1.tex","w");
fprintf(fid12,"Vs V1 0 DC %.11f \n\
R1 V2 V1 %.11f \n\
R2 V3 V2 %.11f \n\
R3 V2 V5 %.11f \n\
R4 0 V5 %.11f \n\
R5 V6 V5 %.11f \n\
R6 V9 V7 %.11f \n\
R7 V7 V8 %.11f \n\
VVc 0 V9 0V \n\
HVc V5 V8 VVc %.11f \n\
GIb V6 V3 V2 V5 %.11f \n\
C1 V6 V8 %.11f", Vs, R1, R2, R3, R4, R5, R6, R7, Kd, Kb, C);


fclose(fid12);

%------  2  -------------------

Vx = V6 - V8;

A2 = [-G1    G1+G2+G3   -G2   -G3   0   0    0    0   0   0   0 ;
        0      -G2       G2    0    0   0    0   -1   0   0   0 ; 
        0       0        0     0    0   G7  -G7   0  -1   0   0 ;
        0       0        0     0    1   0   -1    0   0   0   0 ;
        0      -1        0     1    0   0    0    0   0   1   0 ;
        0       0        0    -1    0   0    1    0   0   0   1 ;
        0       0        0     0    0   0    0    1   0  -Kb  0 ;
        0       0        0     0    0   0    0    0  -Kd  0   1;
        G1+G4 -G1        0    -G4    0   0    0    0   1   0   0;
        1       0        0     0    0   0    0    0   0   0   0 ;
        -G6     0        0     0    0   G6   0    0   1   0   0];
        
B2 = [0;0;0;Vx;0;0;0;0;0;0;0];

NC2=linsolve(A2,B2);

V1_2 = NC2(1);
V2_2 = NC2(2);
V3_2 = NC2(3);
V5_2 = NC2(4);
V6_2 = NC2(5);
V7_2 = NC2(6);
V8_2 = NC2(7);
Ib_2 = NC2(8);
Id_2 = NC2(9);
Vb_2 = NC2(10);
Vd_2 = NC2(11);

Ix = Ib_2 + (V6_2 - V5_2)*G5;
Req = Vx/Ix;
tau = Req * C;


%Criar ficheiro octave teorico
fid21=fopen("teorico2.tex","w+");
fprintf(fid21,"Vx & %f\\\\ \\hline \n\
Ix & %f\\\\ \\hline \n\
Req & %f\\\\ \\hline \n\
Tau & %f\\\\ \\hline ", Vx,Ix,Req,tau)
fclose(fid21)



%ficheiro para ngspice 2
fid22 = fopen("ngspice_2.txt", "w");
fprintf(fid22, "Vs V1 0 DC 0 \n\
R1 V2 V1 %.11f \n\
R2 V3 V2 %.11f \n\
R3 V2 V5 %.11f \n\
R4 0 V5 %.11f \n\
R5 V6 V5 %.11f \n\
R6 V9 V7 %.11f \n\
R7 V7 V8 %.11f \n\
VVd 0 V9 0V \n\
HVd V5 V8 VVd %.11f \n\
GIb V6 V3 V2 V5 %.11f \n\
Vx V6 V8 DC %.11f", R1, R2, R3, R4, R5, R6, R7, Kd, Kb,Vx); 
fclose(fid22);

%------------------------  3  ------------------------



syms v6(t)
t=0:2e-6:20e-3;
v6 = (V8_2 + Vx) * exp(-(t/tau));

Teoria_3 = figure ();
plot(t*1000, v6);
xlabel ("t[ms]");
ylabel ("V_{6n} [V]");
print (Teoria_3, "Teoria_3_Fig.eps","-depsc");
close(Teoria_3);

fid3 = fopen("ngspice_3.txt", "w");
fprintf(fid3, "Vs V1 0 DC 0\n\ 
R1 V2 V1 %.11e \n\
R2 V3 V2 %.11e \n\
R3 V2 V5 %.11e \n\
R4 0 V5 %.11e \n\
R5 V6 V5 %.11e \n\
R6 V9 V7 %.11e \n\
R7 V7 V8 %.11e \n\
VVd 0 V9 0V \n\
HVd V5 V8 VVd %.11e \n\
GIb V6 V3 V2 V5 %.11e \n\
C1 V6 V8 %.11e ic = %.11e \n\
.ic v(V6) = %.11e v(V8) = 0", R1, R2, R3, R4, R5, R6, R7, Kd, Kb,C,Vx, Vx); 
fclose(fid3);

%----------4--------------

f=1000;
Yc = (C*2*pi*f)*i;
Zc = 1/Yc;


A4 = [-G1    G1+G2+G3  -G2   -G3     0    0     0  ;
        0    -G2-Kb     G2    Kb     0    0     0  ; 
        0       Kb       0 -Kb-G5  G5+Yc  0    -Yc ;
        0       0        0    0      0   G6+G7 -G7 ;
        1       0        0    0      0    0     0  ;
        0       0        0    1      0  Kd*G6  -1  ;
        G1     -G1       0   -G4     0  -G6     0  ];

B4 = [0;0;0;0;1;0;0];

NC4=linsolve(A4,B4);

V1_4 = NC4(1);
V2_4 = NC4(2);
V3_4 = NC4(3);
V5_4 = NC4(4);
V6_4 = NC4(5);
V7_4 = NC4(6);
V8_4 = NC4(7);


%Amplitude
AbsV1=abs(V1_4);
AbsV2=abs(V2_4);
AbsV3=abs(V3_4);
AbsV5=abs(V5_4);
AbsV6=abs(V6_4);
AbsV7=abs(V7_4);
AbsV8=abs(V8_4);

%Argumento de cada no
ArgV1=arg(V1_4);
ArgV2=arg(V2_4);
ArgV3=arg(V3_4);
ArgV5=arg(V5_4);
ArgV6=arg(V6_4);
ArgV7=arg(V7_4);
ArgV8=arg(V8_4);

%Faz tabela teorica4 Amplitudes
fid41=fopen("teorico4_Amp.tex","w+");
fprintf(fid41,"V1 & %.11f\\\\ \\hline \n\
V2 & %.11f\\\\ \\hline \n\
V3 & %.11f\\\\ \\hline \n\
V5 & %.11f\\\\ \\hline \n\
V6 & %.11f\\\\ \\hline \n\
V7 & %.11f\\\\ \\hline \n\
V8 & %.11f\\\\ \\hline", AbsV1,AbsV2,AbsV3,AbsV5,AbsV6,AbsV7,AbsV8)
fclose(fid41)


%Faz tabela teorica4 argumentos
fid42=fopen("teorico4_Arg.tex","w+");
fprintf(fid42,"V1 & %.11e\\\\ \\hline \n\
V2 & %.11e\\\\ \\hline \n\
V3 & %.11e\\\\ \\hline \n\
V5 & %.11e\\\\ \\hline \n\
V6 & %.11f\\\\ \\hline \n\
V7 & %.11f\\\\ \\hline \n\
V8 & %.11f\\\\ \\hline",ArgV1,ArgV2,ArgV3,ArgV5,ArgV6,ArgV7,ArgV8)
fclose(fid42)


%----------5------------

t=-5e-3:2e-6:20e-3;

ct=1;

while ct<= length (t) 
  if t(ct)>=0
    Plot_V6(ct)= V6_4*sin(2*pi*f*t(ct)) + Vx*exp(-t(ct)/tau);
    Plot_Vs(ct)= sin(2*pi*f*t(ct));
 elseif t(ct)<0
   Plot_V6(ct) = V6;
   Plot_Vs(ct) = Vs;
 end
 ct=ct+1;
end

 
%figura para Teoria 
Teoria_5 = figure();
plot(t, Plot_V6, t, Plot_Vs);
xlabel ("t");
ylabel ("v_6 [V]    v_s [V]");
legend("v6","vs");
print (Teoria_5, "Teoria_5_Fig.eps","-depsc");
close(Teoria_5);


%ficheiro para ngspice
fid5 = fopen("ngspice_5.txt", "w");
fprintf(fid5, "Vs V1 0 0.0 ac 1.0 sin(0 1 1k) \n\
R1 V2 V1 %.11f \n\
R2 V3 V2 %.11f \n\
R3 V2 V5 %.11f \n\
R4 0 V5 %.11f \n\
R5 V6 V5 %.11f \n\
R6 V9 V7 %.11f \n\
R7 V7 V8 %.11f \n\
VVd 0 V9 0V \n\
HVd V5 V8 VVd %.11f \n\
GIb V6 V3 V2 V5 %.11f \n\
C1 V6 V8 %.11e ic = %.11f \n\
.ic v(V6) = %.11f v(V8) = 0", R1, R2, R3, R4, R5, R6, R7, Kd, Kb,C,Vx, Vx); 
fclose(fid5);




%--------------6-------------
phi_Vs = pi/2
Vs_p= 1*power(e,-j*phi_Vs)


f =-1:0.1:6; %Hz
w = 2*pi*power(10,f);


Vs_p= 1*power(e,-j*phi_Vs);
Zc=1. ./ (j .* w .* C);

B_A = [Kb+1./R2, -1./R2, -Kb, 0;
     1./R3-Kb,  0, Kb-1./R3-1./R4, -1./R6;
     Kb-1./R1-1./R3, 0, 1./R3-Kb, 0;
     0, 0, 1., Kd/R6-R7/R6-1.];
B_B = [0; 0; -Vs_p/R1; 0];

B_C=linsolve(B_A,B_B); % v2, v3, v5, v7
 
V8_6 = R7*(1./R1+1./R6)*B_C(4) + 0*Zc;
V6_6 = ((1./R5+Kb)*B_C(3)-Kb*B_C(1)+ (V8_6 ./ Zc)) ./ (1./R5 + 1. ./ Zc);
Vc = V6_6  - V8_6;
Vs_6 = power(e,j*pi/2) + 0*w;


hf = figure ();
plot (f, 20*log10(abs(Vc)), "r");
hold on;
plot (f, 20*log10(abs(V6_6)), "b");
hold on;
plot (f, 20*log10(abs(Vs_6)), "g");
xlabel ("log_{10}(f) [Hz]");
ylabel ("Magnitude v_s(f), v_6(f), v_c(f) [dB]");
print (hf, "db_teorica.eps", "-depsc");



V6_a = 180/pi*(angle(V6_6));


for  i=1:length(V6_a)
        if(V6_a(i)<=-90) 
		V6_a(i) += 180;
	elseif (V6_a(i)>=90) 
		V6_a(i) -= 180;
endif
endfor


hf = figure ();
plot (f, 180/pi*(angle(Vc) + pi), "r");
hold on;
plot (f, V6_a, "b");
hold on;
plot (f, 180/pi*angle(Vs_6), "g");


xlabel ("log_{10}(f) [Hz]");
ylabel ("Phase v_s(f), v_6(f), v_c(f) [degrees]");
print (hf, "fase_teorica.eps", "-depsc");