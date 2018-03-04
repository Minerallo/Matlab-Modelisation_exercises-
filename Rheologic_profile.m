% programme M2 sciences de la terre
%First program
%simple rheologic profil
clear;
close all
%----------------------------------------------------------------------------
%parameters
T0=0;               %temperature de surface [c]
Tm=1350             %temperature d'asthénosphere [c]
kappa=1e-6;         %diffusivity [m2/s]
zmin=0 ;            %profondeur min  [km]
zmax=200;           %prof max [Km]
dz=1;                 % resolution
Tage=10;
%-----------------------------------------------------------------------------
% unité (SI)
zmin=zmin*1e3                 % en [m]
zmax=zmax*1e3
dz=dz*1e3
Tage=Tage*1e6*365*24*3600;    % (Ma)
%----------------------pour un seul valeur Tage (en 1D)------------------------------------------------------

Tab1d_z=zmin:dz:zmax;                         % remplire un tab par valeur de zmin à zmax avec un pas de dz=1
T= zeros(size(Tab1d_z))                       % faire un vecteur de même size que Tab

T=Tm+(T0-Tm)*erfc(Tab1d_z./(2*sqrt(kappa*Tage)))

T=T+273.15;            % [k]

%---------------------------------
n=3.6;
A=(1e6)^(-n)*10^(4.5);        %Pa/s
Q=535000;                 % j/mol
R=8.3144;           %j/k*mol
eps=1e-12;           % s-1
phi=30;
C=50000;
rho=3300;             % kg/m2
g=9.81;            % m/s
%segma3=rho*g*z
%segma=size(Tab1d_z)
%-----------------------------loi flg------------
a=exp(Q./(R.*T))./A;
tho=zeros(size(Tab1d_z));
tho =(eps.*a).^(1/n);

figure(2),hold on
plot(tho/1e6,-Tab1d_z/(1e3)),
xlim([0 4000]), ylim([-200 0])

p=rho.*g.*Tab1d_z
tho1=p*sin(phi);
figure(2)
plot(tho1/1e6,-Tab1d_z/(1e3)),
xlim([0 4000]), ylim([-200 0])
% %-----------------------------------

% %pour un extension
segma1=rho.*g.*Tab1d_z
segma3=-((2*C*cosd(phi)-segma1.*(1-sind(phi)))/(1+sind(phi)))
tho1=segma1-segma3
hold on
plot(tho1/1e6,-Tab1d_z/(1e3))
ylabel('profondeur [Km]')
xlabel('tho')
