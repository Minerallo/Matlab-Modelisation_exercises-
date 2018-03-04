%geothermie1D_script1
clear ;
close all ;
clc ;

%------------------------------------------------

Ts = 0;          %  Surface Temperature [C]
Tm = 1350;       %  Asth. Temp.         [C]
kappa = 1e-6;    %  Diffusivity         [m²/s]

zmin = 0;         % min Depth            [km]
zmax = 200;       % max Depth            [km]
dz   = 1;         %resolution            [km]

Tage = 50          % Thermal age        [Ma]

%--------------------------------------------------
% Units [Syst, Int.];
zmin = zmin*1e3;    % min Depth             [m]
zmax = zmax*1e3;   % max Depth              [m]
dz   = dz*1e3;
Tage = Tage*1e6*365*24*3600; % Thermal age  [s]
%--------------------------------------------------

Tab1d_z = [zmin:dz:zmax];  % Depths       [km]

T = Tm + (Ts-Tm)*erfc(Tab1d_z./sqrt(kappa*Tage));
figure(1),
plot(T,-Tab1d_z/1e3,'r-','lineWidth',2);
xlabel('Temperature [C]');
ylabel('z [km]');




