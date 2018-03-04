clear all;close all;
%Definition de la grille, condition PT
dP=10^9;dT=100;
P1=0:dP:100*10^9;
T1=100:dT:4000;
[P,T]=meshgrid(P1,T1); 
%Valeur paramètres model
R     = 8.314;          % J/mol, gas constant
Pr    = 100000;         % Pa, reference pressure
Tr    = 298.15;         % K, reference temperature 
Hr    = -601500.00;     % J
Vr    = 1.12228e-5;     % J/Pa, 
B     = 30179500000;    % Pa 
c(1)  = 1.96612;        % dimensionless
c(2)  = 4.12756;        % dimensionless
c(3)  = 0.53690;        % dimensionless
dH(1) = 2966.88;        % J
dH(2) = 5621.69;        % J 
dH(3) = 27787.19;       % J
dV(1) = 3.52971e-8;     % J/Pa
dV(2) = dV(1);          % J/Pa
dV(3) = 1.9849568e-6;   % J/Pa

% Computing Gibbs free energy at given P and T
F=5/4.*(Pr+B)^(1/5).*((P+B).^(4/5)-(Pr+B).^(4/5));
G=Hr+Vr.*F;

for i=1:1:3
    e{i}=exp(-(dH(i)+dV(i).*F)./R./T);
    e0{i}=exp(-dH(i)./R./Tr);
    G=G+c(i).*(R.*T.*log(1-e{i})-dH(i).*e0{i}./(1-e0{i}));
end
%surf(G);colormap(jet);colorbar;%shading interp;
