clear ;
close all ;
clc ;

% Physical parameters-----------------------------------------------------
L = 100;            %   Lenght of modeled domain [m]
Tmagma = 1200;      %   Temperature of magma [C]
Trock = 300; %   Temperature of country rock [C]
W = 50;             %   width of dike [m]
ka = 2;              %  W/m2
kb = 4;
kc = 3;
rhopa = 2700;
rhopb = 3300;
rhopc = 2800;
cpa = 800;
cpb = 1200;
cpc = 900;
%------------------------------------------------------------------------
xmin = -100;
xmax = 100;
zmin = -200;
zmax =200;
nx = 50 ;
nz = 100 ;
dx = (xmax -xmin) / (nx-1) ;
dz = (zmax -zmin) / (nz-1) ;

nt = 100 ;
dt = 100*3600*24 ;
x = [xmin:dx:xmax] ;
z = [zmin:dz:zmax] ;

%-------------------------------------------------------------------------

T = ones(nx,nz).*Trock; % matrice Temperature de taille x de 300ï¿½C

% Bordures limites du dike
xl = 0-W/2;
xr = 0+W/2;
zn = 0-W/2;
zs = 0+W/2;

rhop = zeros (nx,nz)+2700;
cp = zeros (nx,nz)+800;


%Boucle dike bordures
for  i = 1 : nx       % Sur l'ensemble des donnï¿½es
    for j= 1 : nz
        if xl<x(1,i)&& x(1,i)<xr && zn<z(1,j) && z(1,j)<zs            % Si x se trouve dans l'interval du dike
            T(i,j) = Tmagma   ;         % alors x = 1200 ï¿½C
            rhop (i,j) = rhopb ;
            cp (i,j) = cpb ;
        end
        if xl>x(1,i)
            rhop (i,j) = rhopa ;
            cp (i,j) = cpa ;
        end
        if xr<x(1,i)
            rhop (i,j) = rhopc ;
            cp (i,j) = cpc ;
        end
        if z(1,j)<zn
            rhop (i,j) = rhopa;
            cp (i,j) = rhopa;
        end
        if z(1,j)>zs
            rhop (i,j) = rhopc ;
            cp (i,j) = cpc ;
        end
    end
end

t = 0;   % a t=0


k = ones (nx-1,nz-1).*2;
for i=1 : nx-1                            %boucle pour k
    for j=1 : nz-1
        xc(1,i) = (x(1,i)+x(1,i+1))/2;        %nouvelle grid pour k
        zc(1,j) = (z(1,j)+z(1,j+1))/2;
        if xl<xc(1,i) && xc(1,i)<xr && zn<zc(1,j) && zc(1,j)<zs
            k (i,j) = kb ;
        end
    end
end


A = sparse(nx*nz,nx*nz);            % Matrice A

num = 1;                            %Assigne un numéro à chaque valeur
for i=1:nx
    for j=1:nz
        Number(i,j) = num;
        num = num+1;
    end
end
Numt= Number' ;

for i = 2:nx-1             % Difference finie et calcul des coefficients
    for j = 2:nz-1
        sx = dt/(cp(i,j)*rhop(i,j)*(dx*dx));
        sz =dt/(cp(i,j)*rhop(i,j)*(dz*dz));
        A( Number(i,j), Number(i ,j)) = sx*k(i,j)+sx*k(i-1,j)+sz*k(i,j)+sz*k(i,j-1)+1;
        A( Number(i,j), Number(i+1,j )) =-sx*k(i,j);
        A( Number(i,j), Number(i ,j+1)) = -sz*k(i,j);
        A( Number(i,j), Number(i-1,j )) = -sx*k(i-1,j);
        A( Number(i,j), Number(i ,j-1)) = -sz*k(i,j-1);
    end
end

for  j = 1:nz          % 1 au bord de la matrice
    A (Number (1,j) , Number (1, j))= 1 ;
    A (Number (nx,j) , Number (nx, j)) = 1;
end

for  i = 1 : nx
    A (Number(i,1), Number (i,1))= 1 ;
    A (Number(i,nz), Number (i,nz))= 1 ;
end

T1D = zeros (nx*nz,1)
k = 0
for i = 1 : nx
    for j =  1 : nz
        k=k+1 ;
        T1D(k,1) = T(i,j) ;
    end
end


Tnew2D = zeros (nx,nz) ;      %calcul de la nouvelle temperature et boucle sur le temps
for  n=1: nt
    Tnew = A\T1D;
    T1D = Tnew;
    k1=0
    for i = 1 : nx
        for j =  1 :nz
            k1=k1+1;
            Tnew2D (i,j)=   T1D(k1,1);
        end
    end
    
    t =   t+dt;
    
    figure(1), clf
    surf(x,z,Tnew2D');
    shading interp
    %plot(x,Tnew2D,'r-','lineWidth',2);
    ylabel('Temperature [C]');
    xlabel('x [m]');
    title ('Diffusion thermique implicite 2D + k, rho, Cp variables')
    axis equal;
    drawnow
    pause(0.1)
end




