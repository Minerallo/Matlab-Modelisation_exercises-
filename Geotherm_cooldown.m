% programme M2 sciences de la terre
clear; close all;
%==========Les Param�tres======================================================
T0 = 273;                   %temp�rature de surface [k]
Tm = 1850;                  %temp�rature initiale de manteau d�lamin� [k]
rhoc = 2650 ;               %densit� de la cro�te sup [Kg.m-3]
rhoci= 2800                 %densit� de la cro�te inf [Kg.m-3]
rhoml= 3250 ;               %densit� du lithos  [Kg.m-3]
rhoa = 3300                 %densit� du asth�nosph�re    [Kg.m-3]

kc = 2.3;                   % la conductivit� de la cro�te sup [W.m-1K-1]
kci = 2.5;                  % la conductivit� de la cro�te inf[W.m-1K-1]
kml= 3.1;                   % la conductivit� de la lithos [W.m-1K-1]
ka = 3.5;                   % la conductivit�  du asth�nosph�re  [W.m-1K-1]

Cpc    = 800;               %capacit� calorifique de la cro�te sup [J.Kg-1.k-1]
Cpci    = 900;              %capacit� calorifique de la cro�te inf [J.Kg-1.k-1]
cpml   = 1230;              %capacit� calorifique  de la lithos [J.Kg-1.k-1]
cpa    = 1260;              %capacit� calorifique  du asth�nosph�re   [J.Kg-1.k-1]

Rc = 1.7e-6;                %chaleur volumique de la cro�te sup�rieure[W.m-3]
Rci = 0.18e-6;              %chaleur volumique de la cro�te inf[W.m-3]
Rml = 0.02e-6;              %chaleur volumique de la litho[W.m-3]
Ra = 0.005e-6;              %chaleur volumique du asth�nosph�re  [W.m-3]

Qo = 11e-3;                 % Flux de chaleur � la base de l�asth�nosph�re pour g�otherme normal[W.m-2]
Qdo = 3.5e-3                % Flux de chaleur  d� la base de l�asth�nosph�re pour g�otherme d�lamin� [W.m-2]

zc = 0 ;                    %  la profondeur � la surface  [m]
zci  = 15*1e3 ;             %  la profondeur de la cro�te sup�rieur  [m]
zml = 30*1e3 ;              %  la profondeur  la cro�te inf�rieure    [m]
za = 100*1e3 ;              %  la profondeur de manteau lithosph�rique   [m]
zo  = 300*1e3 ;             %  la profondeur de l'asth�nosph�re  [m]

zmax = 300 ;                %profondeur max  [Km]
zmin= 0;                    %Profondeur min  [km]
nz   =  300;                % taille de z

tmin=0;                     %temps min  [s]
tmax=3000;                  %temps max [s]
nt = 1000;

zmin=zmin*1e3;                      % profondeur en [m]
zmax=zmax*1e3;                      % la Profondeur [m]
dz  = (zmax - zmin) / (nz) ;        % r�solution

Ma = 1e6*365*24*3600;
tmin = tmin*1e6*365*24*3600;        % Pour mettre l'unit� de temps en [ Ma ] au lieu de [s]
tmax = 100*tmax*1e6*365*24*3600;    % Pour mettre l'unit� de temps en [ Ma ] au lieu de [s]
dt=(tmax - tmin) / (nt);            % Pas de temps
z= zmin:dz:zmax;                    % le tableau de la profondeur de 0 � 300 km
zinv= z(end:-1:1);
t=tmin:dt:tmax;                     % le tableau de temps[m]

rho = zeros (size (z));             % Tableau pour remplissage de densit�
cp = zeros (size (z));              % Tableau pour remplissage de densit�
k = zeros (size (z));               % Tableau pour remplissage de conductivit� thermique
kappa = zeros (size (z));           % Tableau pour remplissage de diffusivit� thermique
T = zeros(size(z));                 % Tableau pour remplissage de la temp�rature
Td = zeros(size(z));                % Tableau pour remplissage de la temp�rature

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

%===== replissage de diffusivit� =============
for i = 1 : size (z,2)
   kappa (i) = k (i) / (rho (i)*cp(i));
end
%=====Affichage de la densit�,capacit� calorifique et la Conductivit� thermique
figure (1)
subplot(1,3,1);
plot(rho,-z/1e3),
xlabel('Densit�');
ylabel('Z (km)');
title ('Rho [kg.m-3]');

hold on
subplot(1,3,2);
plot(cp,-z/1e3),
xlabel('Capacit� thermique');
ylabel('Z (km)');
title ('Cp [J.kg-1.K-1]');

hold on
subplot(1,3,3);
plot(k,-z/1e3),
xlabel('Conductivit� thermique');
ylabel('Z (km)');
title ('k [W.m-1.K-1]');

%=======================d�termination du g�otherme normale===========
z = zinv ;
zo=300*1e3;                     % � l'asth�nosph�re
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
Q3 = 60e-3;                     % le flux chaleur � la surface : connu
zo=zci

for i= 1:size (z,2)
   if (zci > z(i))
       T(i)=(-Rc/(2*kc))*(z(i)-zo)^2+(Q3/kc)*(z(i)-zo)+T(285)
   end
end

%=======================d�termination du g�otherme d�lamin�===========
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
Qd3 = 60e-3;                                % le flux chaleur � la surface : connu
zo=zci

for i= 1:size (z,2)
   if (zci > z(i))
       Td(i)=(-Rc/(2*kc))*(z(i)-zo)^2+(Qd3/kc)*(z(i)-zo)+Td(285)
   end
end
%===============La boucle de changement du gradient g�othermique avec le temps============

t=0
for i = 1 : nt
   Tnew=zeros(size(T));
   Tnew(1,1)=T(1,1);                                                                    %Conditions aux limites
   Tnew(1,size(T,2))=T(1,size(T,2));                                            %Conditions aux limites
   for j=2:size (z,2)-1;
       Tnew(1,j)=Td(1,j)+(T(1,j)-Td(1,j))*erfc(z(1,i)/(2*sqrt(kappa(1,i)*t)));
   end

   t =   t + dt;

   % La boucle pour d�terminer la profondeur de la lithosph�re � une Tnew donn�e
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
   legend ('Geotherme � Eq','T�d�lamination','T�refroidissement')
   %xlim([0 Tm])
   drawnow
   %pause(0.1)

end
