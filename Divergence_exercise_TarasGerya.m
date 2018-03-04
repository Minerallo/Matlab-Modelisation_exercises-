clear all; close all;
%% paramètres de la boite
W=1*10^6; H=1.5*10^6; % [m] Largeur et longueur du model
nodes=36;
[x,y] = meshgrid(linspace(1,W,nodes),linspace(1,H,nodes));
%% paramètres du modèle 
vxo=10^-9*0.33; %[m/s] valeur initiale des composantes du vecteur vitesse
vyo=10^-9; %[m/s]
%% Equation de la vitesse  
Vx=-vxo*sin(2*pi*x./W).*cos(pi*y./H);
Vy=vyo*cos(2*pi*x./W).*sin(pi*y./H);
%% Calcul des dérivés partielles
Vxdx=-vxo*2*pi/W.*cos(2*pi*x./W).*cos(pi*y./H);
Vydy=vyo*cos(2*pi*x./W).*-sin(pi*y./H)*pi/H;
divV=Vxdx+Vydy;
%% Plot
figure(1);
subplot(2,3,1);pcolor(Vx);colormap(jet);colorbar;shading interp;title('Vx m/s')
subplot(2,3,2);pcolor(Vy);colormap(jet);colorbar;shading interp;title('Vy m/s')
subplot(2,3,3);pcolor(Vxdx);colormap(jet);colorbar;shading interp;title('dvx/dx /s')
subplot(2,3,4);pcolor(Vydy);colormap(jet);colorbar;shading interp;title('dvy/dy m/s')
subplot(2,3,5);pcolor(divV);colormap(jet);colorbar;shading interp;title('div(v) /s')
subplot(2,3,6);quiver(x/1000,y/1000,Vx,Vy);axis([0 W/1000 0 H/1000]);title('champ de vitesse');set(gca,'YDir','reverse');




