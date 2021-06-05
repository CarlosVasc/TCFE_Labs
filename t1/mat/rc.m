close all
clear all


r1=1.0162457409e3;
r2=2.05507529065e3;
r3=3.14574535918e3;
r4=4.03531015178e3;
r5=3.10425245864e3;
r6=2.0191831586e3;
r7=1.03927082306e3;
Va=5.11324505775;
Id=1.00396804542e-3;
Kb=7.27390671914e-3;
Kc=8.11887433018e3;

g1=1/r1;
g2=1/r2;
g3=1/r3;
g4=1/r4;
g5=1/r5;
g6=1/r6;
g7=1/r7;

A=[
  g1+g2+g3   -g2   0     0     0     0      -g3     0;
   -g2-Kb     g2   0     0     0     0       Kb     0;
    Kb        0    g5    0     0     0     -Kb-g5   0;
    0         0    0    -g7  g6+g7  -g6      0      0;
    0         0    0     1   -Kc*g6 Kc*g6   -1      0;
    0         0    0     0    0     1        0      0;
    0         0    0     g7   -g7    0       0      -1;
   -g3        0   -g5    0    0    -g4    g4+g5+g3  1;
  ];

B=[ 0 ; 0; Id ; 0 ; 0 ;-Va; -Id;0];

X=linsolve(A,B);
V1=X(1);
V2=X(2);
V3=X(3);
V4=X(4);
V5=X(5);
V6=X(6);
V7=X(7);
It=X(8);

%fprintf(" V1= %f \n V2= %f \n V3= %f \n V4= %f \n V5= %f \n V6= %f \n V7= %f \n It= %f \n\n",X);


M=[
  r1+r4+r3    -r3        -r4        0      0   ;
     0         1         0         0       -Kb;
  -r4          0      r4+r6+r7-Kc   0      0 ;
  0            0         0        1      0;
  r3          -r3        0          0      1 ;
  ];

R=[ -Va ; 0 ; 0 ; Id ; 0];
%format long g
S=linsolve(M,R);
Ia=S(1);
Ib=S(2);
Ic=S(3);

%fprintf(' Ia= %f \n Ib= %f \n Ic= %f \n Vb= %f \n Vd= %f \n Vg= %f \n',X);
Ir3=Ia-Ib;
Ir5=Id-Ib;
Ir4=Ia-Ic;
Ivc=Id-Ic;


fid=fopen("tabelaNos.tex","w+");
fprintf(fid,"V1 & %f\\\\ \\hline \n\
V2 & %f\\\\ \\hline \n\
V3 & %f\\\\ \\hline \n\
V4 & %f\\\\ \\hline \n\
V5 & %f\\\\ \\hline \n\
V6 & %f\\\\ \\hline \n\
V7 & %f\\\\ \\hline \n\
Ivc & %f\\\\ \\hline", V1,V2,V3,V4,V5,V6,V7,It)

fclose(fid);

fid2=fopen("tabelaMalhas.tex","w+");
fprintf(fid,"Ia & %f\\\\ \\hline \n\
Ib & %f\\\\ \\hline \n\
Ic & %f\\\\ \\hline \n\
Id & %f\\\\ \\hline",Ia,Ib,Ic,Id)

fclose(fid2);


fid3=fopen("tabelaAux.tex","w+");
fprintf(fid,"Ir3 & %f\\\\ \\hline \n\
Ir4 & %f\\\\ \\hline \n\
Ir5 & %f\\\\ \\hline \n\
Ivc & %f\\\\ \\hline",Ir3,Ir4,Ir5,Ivc)

fclose(fid3);
