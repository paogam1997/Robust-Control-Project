%% Parametri
l1=0.5;                                           %[m]
l2=0.4;                                           %[m]
h1=ureal('h1',0.35,'Percentage',5);               %[m] 0.35+-5%
h2=ureal('h2',0.181,'Percentage',10);             %[m] 0.181+-10%
M1=3.25;                                          %[kg]
M2=1.90;                                          %[kg]
I1=0.654;                                         %[kg*m^2]
I2=0.117;                                         %[kg*m^2]
Km1=ureal('Km1',1.08,'Percentage',10);            % 1.08+-10%
Tm1=ureal('Tm1',0.005,'Percentage',20);           %0.005+-20%
Km2=ureal('Km2',0.335,'Percentage',10);           %0.335+-10%
Tm2=ureal('Tm2',0.002,'Percentage',20);           %0.002+-20%
T1=ureal('T1',Tm1.NominalValue/2,'Percentage',5); %Tm1/2+-5%
T2=ureal('T2',Tm2.NominalValue/2,'Percentage',5); %Tm2/2+-5%
c1=6.54*10^(-2);
c2=2.32*10^(-2);
g=9.81;            %m/s^2
dteta0=[0;0];
teta0=[0;0];
s=tf('s');
Gm1_nom=(Km1.NominalValue)/(Tm1.NominalValue*s+1);  %Scelta modelli nominali senza ritardi
Gm2_nom=(Km2.NominalValue)/(Tm2.NominalValue*s+1);
Gm_nom=blkdiag(Gm1_nom,Gm2_nom);  %Rappresentazione a blocchi
%% Matrici Linearizzato
M=[I1+M1*h1^2+M2*l1^2   l1*M2;
    l1*M2               I2+M2*h2^2];
P=[-M1*g  0;
    0     -M2*g];
N=[c1+c2  -c2;
    -c2    c2];
G=[1   0;
  -1   1];
T=[1  -1;
   0   1];
%% Forma di stato modello
A=[zeros(2,2)  eye(2);
    -M\P       -M\N];
B=[zeros(2,2)  zeros(2,2);
   M\G          M\T];    %in=[tau_m;tau_d]
C=[1 0 0 0;
   -1 1 0 0];  %out=[teta1;teta2-teta1]  
D=zeros(2,4);
Plant=minreal(ss(A,B,C,D)); %è un uss
%Un po di plot per inquadrare il problema
bodemag(Gm1_nom,'b',Gm2_nom,'g');grid on; title('Attuatori nominali'); legend('per \tau_m1','per \tau_m2'); figure()
sigma(usample(Plant,50));grid on; title('\sigma di Plant:in=2 coppie, out=angoli relativi');


%% %%%%%%%%%%% LQG senza integratore%%%%%%%%%%
cond_iniz=[teta0(1);teta0(2);dteta0(1);dteta0(2);0;0];
Gnom=Plant.NominalValue; %uniformare notazione
Gmnom=Gm_nom;
%Creo sistema nominale generale
systemnames='Gnom Gmnom'; 
inputvar='[tau_d(2);u(2)]';
outputvar='[Gnom]';
input_to_Gnom='[tau_d;Gmnom]';
input_to_Gmnom='[u]';
sysoutname='P_attuato';
cleanupsysic='yes'; sysic;
P_attuato=minreal(ss(P_attuato)); %Inglobo attuatori nel sistema
%estraggo le matrici dinamiche (!nuovo stato di dim 6 dopo sysic!)
A_an=P_attuato.A;
B_an=P_attuato.B;
C_an=P_attuato.C;
D_an=P_attuato.D;
%Def rumore wd additivo nel modello, non è tau_d
stdv_wd=0.03;
W=blkdiag(stdv_wd^2,stdv_wd^2,stdv_wd^2,stdv_wd^2,stdv_wd^2,stdv_wd^2);
wd=(W^1/2)*randn(size(W,1),1);  %ugualmente wd=out.wd prendendolo da simulink
%Def rumore misura wn
stdv_wn=0.03;
V=blkdiag(stdv_wn^2,stdv_wn^2);
wn=(V^1/2)*randn(2,1);          %ugualmente wn=out.wn prendendolo da simulink
%Risolvo ARE del problema LQR
Q=blkdiag(50*eye(4),eye(2)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%da modificare, peso su stati
R=10*eye(2); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%da modificare, peso su ingressi
[X,Kr,cl_eig]=icare(A_an,B_an(:,3:4),Q,R,[],[],[]); %X=soluzione ARE, Kr=retoazione ottima u=-Kr*x
%Risolvo ARE filtro di Kalman
[Ytilde,~,cl_eig1]=icare(A_an',C_an',W,V,[],[],[]); %risultati parziali, vanno sistemati per tornare nella forma del mio problema
Y=Ytilde';
Kf=Y*C_an'/V;
K_LQG=[A_an-B_an(:,3:4)*Kr-Kf*C_an    Kf;
                   -Kr                zeros(size(Kr,1),size(Kf,2))] %forma brutta per forma di stato
%% %%%%%%%%%%%%%%%%% LQG con integratore %%%%%%%%%%%%%%%%%%%%%%% %%
cond_iniz=[teta0(1);teta0(2);dteta0(1);dteta0(2);0;0];
Gnom=Plant.NominalValue; %uniformare notazione
Gmnom=Gm_nom;
%Creo sistema nominale generale
systemnames='Gnom Gmnom';
inputvar='[tau_d(2);u(2)]';
outputvar='[Gnom]';
input_to_Gnom='[tau_d;Gmnom]';
input_to_Gmnom='[u]';
sysoutname='P_attuato';
cleanupsysic='yes'; sysic;
P_attuato=minreal(ss(P_attuato)); %Inglobo attuatori nel sistema
%estraggo le matrici dinamiche (!nuovo stato di dim 6 dopo sysic!)
A_an=P_attuato.A;
B_an=P_attuato.B;
C_an=P_attuato.C;
D_an=P_attuato.D;
%Def disturbo wd=B_an(:,1:2)*tau_d
stdv_wd=0.03;
W=blkdiag(stdv_wd^2,stdv_wd^2,stdv_wd^2,stdv_wd^2,stdv_wd^2,stdv_wd^2);
wd=(W^1/2)*randn(size(W,1),1);  %ugualmente wd=out.wd prendendolo da simulink
%Def rumore misura wn
stdv_wn=0.03;
V=blkdiag(stdv_wn^2,stdv_wn^2);
wn=(V^1/2)*randn(2,1);          %ugualmente wn=out.wn prendendolo da simulink
[n,m]=size(B_an(:,3:4));     %n=#stati,m=#in
%Genero l'impianto aumentato
Aaum=[A_an   zeros(n,m);
      -C_an  zeros(m,m)];
Baum=[B_an(:,3:4);
       -D_an(:,3:4)];
%Definisco i pesi
Q=[50*eye(n,n) zeros(n,m);zeros(m,n) 50*eye(m,m)]; %peso su errore integrato
R=10*eye(m); %peso su ingresso
Kr=lqr(Aaum,Baum,Q,R); %retroaz ottima lqr
Krp=Kr(1:m,1:n); Kri=Kr(1:m,n+1:n+m); %separo parte relativa a stato dell'integratore (questo non devo stimarlo)
Bnoise=eye(n); %modello di wd, è additivo
Kf=lqe(A_an,Bnoise,C_an,W,V); %guadagno del filtro di Kalman
%Reinserisco lo stato di integratore
Ac=[zeros(m,m)            zeros(m,n);
    -B_an(:,3:4)*Kri      A_an-B_an(:,3:4)*Krp-Kf*C_an];
Bc=[eye(m)     -eye(m);
    zeros(n,m)  Kf];
Cc=[-Kri -Krp];  Dc=[zeros(m,m) zeros(m,m)];
Klqg_int=ss(Ac,Bc,Cc,Dc); %Controllore finale lqg con integratore


%% %%%%%%%%%%% Modellazione Incertezze %%%%%%%%%%%%%%%
%% Attuatori 
% Modello incertezze moltiplicative in ingresso
omega=logspace(-2,6,100);
for k=0.972:0.01:1.188  %Primo attuatore
    for Tm=1/250 :0.01: 3/500
        for t=19/800 :0.01: 21/800
            Gm1_p= (k*exp(-t*s))/(Tm*s+1);
            rel_error1=(frd(Gm1_p,omega)-frd(Gm1_nom,omega))/(frd(Gm1_nom,omega));
            bodemag(rel_error1,'k--',omega);grid on
            hold on
            title('Errore relativo attuatore 1')
        end
    end
end
%% Trovo l'upper bound Wm1
[freq,resp_db]=ginput(10);
resp=10.^(resp_db/20);
sys=frd(resp,freq);
Wm1=fitmagfrd(sys,2,0); %ordine 2, no ord rel
Wm1=tf(Wm1)  %Wm1=tf([2.44 102.6 387.9],[1 109.9 3795])
hold on
bode(Wm1);
hold off
%% Secondo attuatore
for k=0.3015:0.01:0.3685  
    for Tm=1/625 :0.01: 3/1250
        for t=19/20000 :0.01: 21/20000
            Gm2_p= (k*exp(-t*s))/(Tm*s+1);
            rel_error2=(frd(Gm2_p,omega)-frd(Gm2_nom,omega))/(frd(Gm2_nom,omega));
            bodemag(rel_error2,'k--',omega);grid on
            hold on
            title('Errore relativo attuatore 2')
        end
    end
end
%% Trovo l'upper bound Wm2
[freq,resp_db]=ginput(10);
resp=10.^(resp_db/20);
sys=frd(resp,freq);
Wm2=fitmagfrd(sys,2,0); %ordine 2, no ord rel
Wm2=tf(Wm2)    %Wm2=tf([2.632 1867 2.081*10^5],[1 2260 2*10^6]) o    Wm2=tf([2.61 1349 1.522*10^5],[1 2408 1.528*10^6])
hold onbode(Wm2);
hold off
%% Attuatori reali
Delta_att1=ultidyn('Delta_att1',[1 1]);
Gm1=Gm1_nom*(1+Wm1*Delta_att1);
Delta_att2=ultidyn('Delta_att2',[1 1]);
Gm2=Gm2_nom*(1+Wm2*Delta_att2);
Wm=blkdiag(Wm1,Wm2);
Delta_att=blkdiag(Delta_att1,Delta_att2);
Gatt=Gm_nom*(eye(2)+Wm*Delta_att);
figure();
bodemag(usample(Gm1,50),'b',usample(Gm2,50),'g');grid on ;title('attuatori reali'); legend('G_m1','G_m2');
%% Creo sistema uss generale y=Gsys*[tau_d;u]
systemnames='Plant Gatt';
inputvar='[tau_d(2);u(2)]'; %rumore di misura è messo dopo! Questo è il sistema con attuatori
outputvar='[-Plant]'; 
input_to_Plant='[tau_d;Gatt]';
input_to_Gatt='[u]';
sysoutname='Gsys';
cleanupsysic='yes';sysic;


%% Def i pesi per prestazioni
%Coloro il rumore bianco di misura, lo concentro  in AF
wn=(2*10^-5)*(10*s+1)/(0.1*s+1);
Wn=blkdiag(wn,wn); 
%Peso su S per il tracking
Mi=2;  wbi=1*10^1;   Ai=10^-6;   enne=1;
wp= (s/Mi^(1/enne) +wbi)^enne/(s+wbi*Ai^(1/enne))^enne;
Wp=blkdiag(wp,wp);
%Wu passa alto per ridurre modulo di ingresso , peso di KS
wu=6*10^-5*(s)/(s/10^1 +1);
Wu=blkdiag(wu,wu);
%Wt passa alto per ridurre sensitività a rumore misura, peso di T
wt=(0.5*10^0)*(s/0.1+1)/(s/10+1);
Wt=blkdiag(wt,wt);
% bodemag(1/wp,'b',1/wt,'g'); grid on; legend('1/wp','1/wt'); title('funzioni peso');
bodemag(1/wp,'b',1/wu,'r',1/wt,'g'); grid on; legend('1/wp','1/wu','1/wt'); title('funzioni peso');
%% Creo sistema generalizzato P in forma di Doyle con i pesi interni
systemnames='Gsys Wp Wu Wt  Wn';  %
inputvar='[ni(2);tau_d(2);u(2)]';
outputvar='[Wp;Wu;Wt;-Wn-Gsys]';   %
input_to_Gsys='[tau_d;u]';
input_to_Wp='[Wn+Gsys]';
input_to_Wt='[Gsys]';
input_to_Wu='[u]';
input_to_Wn='[ni]';
sysoutname='Pdoyle';
cleanupsysic='yes'; sysic;


%% %%%%%%%%%%%%%%%%%%%%%% Sintesi Hinf %%%%%%%%%%%%%%%%5
cond_iniz=[teta0(1);teta0(2);dteta0(1);dteta0(2)];
%Estraggo valore nominale Pnom usato per sintesi Hinfinito
Pnom=Pdoyle.NominalValue;
nmeas = 2; %dim ingresso controllore
ncon = 2; %dim uscita controllore
gmin = 0; %minimo di gamma 
gmax = 10; %max di gamma
opts = hinfsynOptions('Display','on');
%Genero controllore
[K_hin,clp,gamma,info_hinf] = hinfsyn(Pnom,nmeas,ncon,[gmin,gmax],opts);
if gamma>=1
    disp('non rispetta pienamente le specifiche');
else
    disp('Rispetta le specifiche')
end
sigma(clp,'b',ss(gamma),'r'), grid on;
title('Valori singolari clp')
K = K_hin;
Ktf=minreal(zpk(tf(K)));
%Controllare risposte a cl
%risposte ai disturbi
%% Studio Robustezza con mu-analisi
% Definisco CL incerto
clp_ic = lft(Pdoyle,K);
%Ne prendo il ufrd
omega = logspace(-2,5,100);
clp_g = ufrd(clp_ic,omega);
%Calcolo autovaori a CL così da verificare la Nominale Stabilità del
%sistema
format long;
autovalori_CL=eig(clp_ic.NominalValue) 
[i,j,v]=find(autovalori_CL>=0); %Queste 7 righe fanno un check automatico della presenza di autovalori instabili
 if isempty(v) 
disp('OK autovalori tutti negativi, hai Nominale Stabilità')
else
disp('FERMO, hai un autovalore instabile in v')
v
end

% Analisi di Robusta stabilità
opt = robopt('Display','on');
[stabmarg,destabu,report_RS,info] = robuststab(clp_g,opt);
report_RS
figure(1)
loglog(info.MussvBnds(1,1),'r-',info.MussvBnds(1,2),'b--')
grid
title('Robust stability')
xlabel('Frequency (rad/s)')
ylabel('\mu')
legend('\mu-upper bound','\mu-lower bound')

% Prestazioni Nominali
figure(2)
sv = sigma(clp_ic.Nominal,omega);
sys_frd = frd(sv(1,:),omega);
bodemag(sys_frd,'r-')
grid
title('Nominal performance')
xlabel('Frequency (rad/s)')

% Prestazioni Robuste
opt = robopt('Display','on')
[perfmarg,perfmargunc,report_RP,info] = robustperf(clp_g,opt);
report_RP
figure(3)
semilogx(info.MussvBnds(1,1),'r-',info.MussvBnds(1,2),'b--')
grid
title('Robust performance')
xlabel('Frequency (rad/s)')
ylabel('\mu')
legend('\mu-upper bound','\mu-lower bound')


%% %%%%%% MU-SINTESI, D-K ITERATION %%%%%%%%%%%
cond_iniz=[teta0(1);teta0(2);dteta0(1);dteta0(2)];
omega = logspace(-3,5,100);
nmeas = 2; %dim ingresso controllore
ncont = 2; %dim uscita controllore

opts=musynOptions('Display','full','MaxIter',100,'TolPerf',0.001,'FrequencyGrid',omega)
[K_mu,CLPperf,info_mu]=musyn(Pdoyle,nmeas,ncont,opts);
K=K_mu;  %step 4 ok, il 20 è una goduria


%% K_4 DK 

k11= tf([],[]);
k12= tf([],[]);
k21= tf([],[]);
k22= tf([],[]);
K=[k11 k12;k21 k22];

%% K_20  DK

k11= tf([],[]);
k12= tf([],[]);
k21= tf([],[]);
k22= tf([],[]);
K=[k11 k12;k21 k22];