% poisson 2D
clear all; close all;

xsize=1000; %km
ysize=1500;

xpnt=31;
ypnt=41;

dx=xsize/(xpnt-1);
dy=ysize/(ypnt-1);

%sparse matrix de coefficient
L=sparse(xpnt*ypnt,xpnt*ypnt);
R=zeros(xpnt*ypnt,1);

for i=1:1:ypnt
    for j=1:1:xpnt
        k=(j-1)*ypnt+i %global index de l'inconnu
        %pour les bordure
        if (i==1 || i==ypnt || j==1 || j==xpnt)
            L(k,k)=1;
            R(k,1)=0;
        else
            L(k,k-ypnt)=1/dx^2;
            L(k,k+xpnt)=1/dx^2;
            L(k,k)=-2/dx^2-2/dy^2;
            L(k,k+1)=1/dy^2;
            L(k,k-1)=1/dy^2;
            R(k,1)=1;
        end
    end
end

%vecteur solution S
S=L\R;

% make(2D)
FI =zeros(ypnt,xpnt);
for i=1:1:ypnt
    for j=1:1:xpnt
        k=(j-1)*ypnt+i;
        FI(i,j)=S(k);
    end
end

x=0:dx:xsize;
y=0:dy:ysize;

figure(1);
surf(x,y,FI);
light;
shading interp;
colorbar;
lighting phong;
title('Solution of 2D Poisson equation')
xlabel('x, km')
ylabel('y, km')
zlabel('Gravity potential, J/kg')


