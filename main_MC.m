clear
clc
format long g

% parametri modela
p.S0=100;
p.K = 100;
p.T=1;
p.r = 0.02;

p.f0=0.15;
p.f1 = -0.1;
p.x0=0;

% casovni koraki
Nlist=[2 4 6 8 10 12 16];

% Monte Carlo
M=1000000;
batchSize = 10000;

% referencna vrednost
N_ref=1024;
M_ref = 1000000;

% izracun referencne vrednosti
refFile=sprintf( ...
 'reference_extrapEuler2_S0%d_K%d_r%.4f_f0%.4f_f1%.4f_N%d_M%d.mat', ...
 p.S0,p.K,p.r,p.f0,p.f1,N_ref,M_ref);

if isfile(refFile)
 load(refFile,'ref')
else
 rng(1001,'twister')
 ref=mc_reference_extrap_euler2(p,N_ref,M_ref,batchSize);
 save(refFile,'ref')
end


% tabela referencne vrednosti
Table3=table(ref.N_ref,ref.delta_ref,ref.M_ref,ref.V_ref,ref.SE_ref, ...
 'VariableNames',{'N_ref','delta_ref','M_ref','V_ref','SE_ref'});

disp(' ')
disp('Referencna vrednost')
disp(Table3)

% obe metodi
nN=numel(Nlist);
delta=p.T./Nlist(:);

V_ST1_ant=zeros(nN,1);
SE_ST1_ant = zeros(nN,1);
SE_ST1_brez=zeros(nN,1);

SEdiff_ST1=zeros(nN,1);
E_ST1 = zeros(nN,1);
Ratio_ST1=zeros(nN,1);

V_SW2_ant=zeros(nN,1);
SE_SW2_ant = zeros(nN,1);
SE_SW2_brez=zeros(nN,1);

SEdiff_SW2=zeros(nN,1);
E_SW2 = zeros(nN,1);
Ratio_SW2=zeros(nN,1);

rng(1)

for i=1:nN

 N=Nlist(i);

 fprintf('N = %d\n',N)

 % poenostavljena sibka Eulerjeva metoda
 st1=mc_simplified_euler1(p,N,M,batchSize);

 V_ST1_ant(i)=st1.V_ant;

 SE_ST1_ant(i)=st1.SE_ant;
 SE_ST1_brez(i) = st1.SE_brez;

 E_ST1(i)=abs(st1.V_ant-ref.V_ref);

 SEdiff_ST1(i)=sqrt(st1.SE_ant^2+ref.SE_ref^2);
 Ratio_ST1(i)=E_ST1(i)/SEdiff_ST1(i);


 % poenostavljena sibka Taylorjeva metoda reda 2
 sw2=mc_simplified_weak2(p,N,M,batchSize);

 V_SW2_ant(i)=sw2.V_ant;

 SE_SW2_ant(i)=sw2.SE_ant;
 SE_SW2_brez(i) = sw2.SE_brez;

 E_SW2(i)=abs(sw2.V_ant-ref.V_ref);

 SEdiff_SW2(i)=sqrt(sw2.SE_ant^2+ref.SE_ref^2);
 Ratio_SW2(i)=E_SW2(i)/SEdiff_SW2(i);

end


%% glavna primerjava

Table4=table(Nlist(:),delta,V_ST1_ant,E_ST1,Ratio_ST1, ...
 V_SW2_ant,E_SW2,Ratio_SW2, ...
 'VariableNames',{'N','delta','V_ST1_ant','E_ST1', ...
 'E_ST1_over_SEdiff','V_SW2_ant','E_SW2','E_SW2_over_SEdiff'});

disp(' ')
disp('Primerjava metod')
disp(Table4)


%% primerjava standardnih napak

Table5=table(Nlist(:),delta,SE_ST1_brez,SE_ST1_ant, ...
 SE_ST1_brez./SE_ST1_ant,SE_SW2_brez,SE_SW2_ant, ...
 SE_SW2_brez./SE_SW2_ant, ...
 'VariableNames',{'N','delta','SE_ST1_brez','SE_ST1_ant', ...
 'SE_ST1_brez_over_SE_ant','SE_SW2_brez','SE_SW2_ant', ...
 'SE_SW2_brez_over_SE_ant'});

disp(' ')
disp('Primerjava standardnih napak')
disp(Table5)


%% empiricni red konvergence

[orderST1,usedST1]=empirical_order(delta,E_ST1,Ratio_ST1,Nlist);
[orderSW2,usedSW2]=empirical_order(delta,E_SW2,Ratio_SW2,Nlist);

Metoda={ ...
 'Poenostavljena sibka Eulerjeva metoda'; ...
 'Poenostavljena sibka Taylorjeva metoda reda 2'};

Teoreticni_sibki_red=[1;2];
Ocenjeni_red=[orderST1;orderSW2];

Uporabljene_vrednosti_N={usedST1;usedSW2};

Table8=table(Metoda,Teoreticni_sibki_red,Ocenjeni_red, ...
 Uporabljene_vrednosti_N);

disp(' ')
disp('Empiricni red konvergence')
disp(Table8)


%% shranjevanje rezultatov

save('localvol_simplified_methods_results.mat', ...
 'p','Nlist','M','batchSize','N_ref','M_ref','ref', ...
 'Table3','Table4','Table5','Table8');


%% lokalna funkcija

function [ord,used]=empirical_order(delta,E_delta,ratio,Nlist)

% samo dovolj zanesljive tocke
idx=(ratio>4) & (E_delta>0) & isfinite(E_delta) & isfinite(ratio);

if nnz(idx)>=2

 c=polyfit(log(delta(idx)),log(E_delta(idx)),1);
 ord=c(1);

 used=strjoin(arrayfun(@num2str,Nlist(idx), ...
     'UniformOutput',false),', ');
else
 ord=NaN;
 used='premalo zanesljivih tock';
end

end