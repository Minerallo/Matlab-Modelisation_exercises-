%geothermie1D_script1
clear ;
close all ;
clc ;

%Physical parameter------------------------------------------------

Ts = 0;          %  Surface Temperature [C]
Tm = 1350;       %  Asth. Temp.         [C]
kappa = 1e-6;    %  Diffusivity         [m²/s]

zmin = 0;         % min Depth            [km]
zmax = 200;       % max Depth            [km]
dz   = 1;         %resolution            [km]

tmin = 0             % min Thermal age        [Ma]
tmax = 50            % max Thermal age        [Ma]
dt = 1             % resolution

%--------------------------------------------------
% Units [Syst, Int.];
zmin = zmin*1e3;   % min Depth              [m]
zmax = zmax*1e3;   % max Depth              [m]
dz   = dz*1e3;
tmin = tmin*1e6*365*24*3600; % minThermal age  [s]
tmax = tmax*1e6*365*24*3600; % maxThermal age  [s]
dt   = dt*1e6*365*24*3600;
%--------------------------------------------------
Tab1d_t = tmin:dt:tmax;  % Temperature  [C]
Tab1d_z = zmin:dz:zmax;  % Depths       [km]

T = zeros(size(Tab1d_z,2),size(Tab1d_t,2));

%Boucle T en fonction de z
for i=1 : size(Tab1d_z,2);
    for j = 1 : size(Tab1d_t,2);
        
        T(i,j) = Tm + (Ts-Tm)*erfc(Tab1d_z(1,i)./sqrt(kappa*Tab1d_t(1,j)));
    end
end
figure(1),
pcolor(Tab1d_t/(1e6*365*24*3600),-Tab1d_z/1e3,T)
shading interp
xlabel('Temps [Ma]');
ylabel('z [km]');
title( 'Geotherme 1D + T')



