% programme M2 sciences de la terre
clear; close all;
%==========Les Paramètres======================================================
T0 = 273;                   %température de surface [k]
Tm = 1850;                  %température initiale de manteau délaminé [k]
rhoc = 2650 ;               %densité de la croûte sup [Kg.m-3]
rhoci= 2800                 %densité de la croûte inf [Kg.m-3]
rhoml= 3250 ;               %densité du lithos  [Kg.m-3]
rhoa = 3300                 %densité du asthénosphère    [Kg.m-3]

kc = 2.3;                   % la conductivité de la croûte sup [W.m-1K-1]
kci = 2.5;                  % la conductivité de la croûte inf[W.m-1K-1]
kml= 3.1;                   % la conductivité de la lithos [W.m-1K-1]
ka = 3.5;                   % la conductivité  du asthénosphère  [W.m-1K-1]

Cpc    = 800;               %capacité calorifique de la croûte sup [J.Kg-1.k-1]
Cpci    = 900;              %capacité calorifique de la croûte inf [J.Kg-1.k-1]
cpml   = 1230;              %capacité calorifique  de la lithos [J.Kg-1.k-1]
cpa    = 1260;              %capacité calorifique  du asthénosphère   [J.Kg-1.k-1]

Rc = 1.7e-6;                %chaleur volumique de la croûte supérieure[W.m-3]
Rci = 0.18e-6;              %chaleur volumique de la croûte inf[W.m-3]
Rml = 0.02e-6;              %chaleur volumique de la litho[W.m-3]
Ra = 0.005e-6;              %chaleur volumique du asthénosphère  [W.m-3]

Qo = 11e-3;                 % Flux de chaleur à la base de l’asthénosphère pour géotherme normal[W.m-2]
Qdo = 3.5e-3                % Flux de chaleur  dà la base de l’asthénosphère pour géotherme délaminé [W.m-2]

zc = 0 ;                    %  la profondeur à la surface  [m]
zci  = 15*1e3 ;             %  la profondeur de la croûte supérieur  [m]
zml = 30*1e3 ;              %  la profondeur  la croûte inférieure    [m]
za = 100*1e3 ;              %  la profondeur de manteau lithosphérique   [m]
zo  = 300*1e3 ;             %  la profondeur de l'asthénosphère  [m]

zmax = 300 ;                %profondeur max  [Km]
zmin= 0;                    %Profondeur min  [km]
nz   =  300;                % taille de z

tmin=0;                     %temps min  [s]
tmax=3000;                  %temps max [s]
nt = 1000;

zmin=zmin*1e3;                      % profondeur en [m]
zmax=zmax*1e3;                      % la Profondeur [m]
dz  = (zmax - zmin) / (nz) ;        % résolution

Ma = 1e6*365*24*3600;
tmin = tmin*1e6*365*24*3600;        % Pour mettre l'unité de temps en [ Ma ] au lieu de [s]
tmax = 100*tmax*1e6*365*24*3600;    % Pour mettre l'unité de temps en [ Ma ] au lieu de [s]
dt=(tmax - tmin) / (nt);            % Pas de temps
z= zmin:dz:zmax;                    % le tableau de la profondeur de 0 à 300 km   
zinv= z(end:-1:1);
t=tmin:dt:tmax;                     % le tableau de temps[m]   

rho = zeros (size (z));             % Tableau pour remplissage de densité
cp = zeros (size (z));              % Tableau pour remplissage de densité
k = zeros (size (z));               % Tableau pour remplissage de conductivité thermique 
kappa = zeros (size (z));           % Tableau pour remplissage de diffusivité thermique 
T = zeros(size(z));                 % Tableau pour remplissage de la température 
Td = zeros(size(z));                % Tableau pour remplissage de la température 

%=============remplissage des tableaux de  rho et cp et k======
for i=1:size(z,2) 
   if (zci > z(i))                 % entre 0 et 15 Km
       rho (i) = rhoc ;
       cp (i)  = Cpc ;
       k (i)   = kc ;
   end
   if (zci <= z(i)) && (zml > z(i))% entre 15 et 30 Km
       rho (i) = rhoci ;
       cp (i)  = Cpci ;
       k (i)   = kci ;
   end
   if (zml<= z(i))&& (za >= z(i))  % entre 30 et 100 Km
       rho (i) = rhoml ;
       cp (i)  = cpml ;
       k (i)   = kml ;
   end
   if (z(i)> za)                   % entre 100 et 300 Km
       rho (i) = rhoa ;
       cp (i)  = cpa ;
       k (i)   = ka ;
   end
end

%===== replissage de diffusivité =============
for i = 1 : size (z,2)
   kappa (i) = k (i) / (rho (i)*cp(i));
end
%=====Affichage de la densité,capacité calorifique et la Conductivité thermique
figure (1)
subplot(1,3,1); 
plot(rho,-z/1e3),
xlabel('Densité');
ylabel('Z (km)');
title ('Rho [kg.m-3]');

hold on
subplot(1,3,2); 
plot(cp,-z/1e3),
xlabel('Capacité thermique');
ylabel('Z (km)'); 
title ('Cp [J.kg-1.K-1]');

hold on
subplot(1,3,3); 
plot(k,-z/1e3),
xlabel('Conductivité thermique');
ylabel('Z (km)'); 
title ('k [W.m-1.K-1]');

%=======================détermination du géotherme normale=========== 
z = zinv ;
zo=300*1e3;                     % à l'asthénosphère
for i=1:size(z,2)
   if (z(i)> za)          
       T(i)=(-Ra/(2*ka))*(z(i)-zo)^2+(Qo/ka)*(z(i)-zo)+Tm
   end
end

Q1=-(-ka*(T(199)-T(200))/1000)    % le flux chaleur entre la profondeur 199 Km et 200 Km
zo=za

for i= 1:size (z,2)
   if (zml<= z(i))&& (za >= z(i))
       T(i)=(-Rml/(2*kml))*(z(i)-zo)^2+(Q1/kml)*(z(i)-zo)+T(199)
   end
end

Q2=-(-kml*(T(270)-T(271))/1000)  % le flux chaleur entre la profondeur 270 Km et 271 Km 
zo=zml

for i= 1:size (z,2)
   if (zci <= z(i)) && (zml > z(i))
       T(i)=(-Rci/(2*kci))*(z(i)-zo)^2+(Q2/kci)*(z(i)-zo)+T(270)
   end
end

%Q3=-(-kc*(T(285)-T(286))/1000)
Q3 = 60e-3;                     % le flux chaleur à la surface : connu
zo=zci

for i= 1:size (z,2)
   if (zci > z(i))
       T(i)=(-Rc/(2*kc))*(z(i)-zo)^2+(Q3/kc)*(z(i)-zo)+T(285)
   end
end

%=======================détermination du géotherme délaminé=========== 
z = zinv ;
zo=300*1e3;
ka=1000
for i=1:size(z,2)
   if (z(i)> za)
       Td(i)=(-Ra/(2*ka))*(z(i)-zo)^2+(Qdo/ka)*(z(i)-zo)+Tm
   end
end

Qd1=-(-ka*(Td(199)-Td(200))/1000)
zo=za

for i= 1:size (z,2)
   if (zml<= z(i))&& (za >= z(i))
       Td(i)=(-Rml/(2*kml))*(z(i)-zo)^2+(Qd1/kml)*(z(i)-zo)+Td(199)
   end
end

Qd2=-(-kml*(Td(270)-Td(271))/1000)
zo=zml

for i= 1:size (z,2)
   if (zci <= z(i)) && (zml > z(i))
       Td(i)=(-Rci/(2*kci))*(z(i)-zo)^2+(Qd2/kci)*(z(i)-zo)+Td(270)
   end
end

%Qd3=-(-kc*(Td(285)-Td(286))/1000)
Qd3 = 60e-3;                                % le flux chaleur à la surface : connu
zo=zci

for i= 1:size (z,2)
   if (zci > z(i))
       Td(i)=(-Rc/(2*kc))*(z(i)-zo)^2+(Qd3/kc)*(z(i)-zo)+Td(285)
   end
end
%===============La boucle de changement du gradient géothermique avec le temps============

t=0
for i = 1 : nt
   Tnew=zeros(size(T));
   Tnew(1,1)=T(1,1);                                                                    %Conditions aux limites
   Tnew(1,size(T,2))=T(1,size(T,2));                                            %Conditions aux limites
   for j=2:size (z,2)-1;
       Tnew(1,j)=Td(1,j)+(T(1,j)-Td(1,j))*erfc(z(1,i)/(2*sqrt(kappa(1,i)*t)));
   end
   
   t =   t + dt;
   
   % La boucle pour déterminer la profondeur de la lithosphère à une Tnew donnée 
   tlitho = 1570 ;
   dif   = abs(Tnew-tlitho);
   match  = dif == min(dif);
   idx = find(dif == min(dif));
   zlitho = z(idx);
   
   figure(2)
   plot(T,-z/1e3,'g','lineWidth',2)
   hold on
   figure(2)
   plot (Td,-z/1e3,'r','lineWidth',2)
   hold on
   ylabel('z [km]');
   xlabel('Temperature [k]');
   drawnow
   figure(2)
   plot(Tnew,-z/1e3,'b','lineWidth',2)
   hold on
   ylabel('z [km]');
   xlabel('Temperature [k]');
   title([' Gradient geothermique, time : ',num2str(t/Ma),' Ma'])
   legend ('Geotherme à Eq','T°délamination','T°refroidissement')
   %xlim([0 Tm])
   drawnow
   %pause(0.1)

end
