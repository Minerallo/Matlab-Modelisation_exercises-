clear ;
close all ; 
clc ; 

% Physical parameters----------------------------------------------------- 
 L = 100;            %   Lenght of modeled domain [m]
 Tmagma = 1200;      %   Temperature of magma [C]
 Trock = 300;        %   Temperature of country rock [C]
 kappa = 1e-6;       %   Thermal doffusivity of rock [m²/s]
 W = 50;             %   width of dike [m]
 %------------------------------------------------------------------------
 xmin = -100;
 xmax = 100;
 dx = 2 ;
 nstep = 100 ;
 dt = 10*3600*24
 x = [xmin:dx:xmax] ;
  %-------------------------------------------------------------------------
% T = zeros(size(x))
%  for i = 1 : size(x,2)
%  T(1,i) = Trock;
%  end

T = ones(size(x)).*Trock; % matrice Temperature de taille x de 300C

% Bordures limites du dike
xl = 0-W/2;
xr = 0+W/2;

%Boucle dike bordures
 for  i = 1 : size (x,2);       % Sur l'ensemble des donnees
    if xl<x(1,i) && x(1,i)<xr ; % Si x se trouve dans l'interval du dike
     T(1,i) = Tmagma            % alors x = 1200 C
 end
 end

% figure(1), clf
% plot(x,T,'r.','lineWidth',2);
% ylabel('Temperature [C]');
% xlabel('x [m]');
 
 
t = 0                           %a t=0
                      % appelons T ou Told  correspond au profil de base

%Boucle de T en fonction du temps
for n=1:nstep                             % Pour chaque interval de temps
    Tnew = zeros (1,size (T,2));          % Je recalcul une nouvelle matrice Tnew
    Tnew (1,1) = T(1,1);              % dont la premiere valeur est egale a Told
    Tnew(1,size(T,2)) = T(1,size(T,2)); % temperature de depart definie au frontiere
        for i=2 : size (x,2)-1            % entre la deuxieme et avant derniere valeur
        cfl = kappa*dt/(dx*dx);
        Tnew (i) = T(i)+ cfl*(T(i+1)-2*T(i)+T(i-1));        %les zeros de la matrice son remplis 
        end

T  =   Tnew;    
t =   t+dt;

figure(1), clf
plot(x,Tnew,'r-','lineWidth',2);
ylabel('Temperature [C]');
xlabel('x [m]');
title ('Diffusion thermique explicite 1D')
ylim([0 Tmagma])
drawnow
pause(0.1)
end

%-------------------------
    