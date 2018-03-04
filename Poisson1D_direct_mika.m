clear all; close all; 
%exos Poisson1D

distance=1e+6; %distance en m du premier point au dernier
npt=1000;%nombre de point 
L=sparse(npt,npt); %est le coefficient devant chaque FI 
%(FI(i,i-1)+FI(i,i+1)-2FI(i,i))/dx^2
dx=distance/(npt-1);
R=zeros(npt,1);

%premier point
L(1,1)=1;
R(1,1)=0;

%point intermédiaire
for i=2:1:npt-1
    %definition des coefficients L
    L(i,i-1)=1/dx^2;
    L(i,i+1)=1/dx^2;
    L(i,i)=-2/dx^2;
    R(i,1)=1;
end

%dernier point
L(npt,npt)=1;
R(npt,1)=0;
%solution vecteur FI
FI=L\R;

%pour plot
x=0:dx:distance;

% Open figure
figure(1);
% Plot the solution
plot(x/1000,FI);
xlabel('x, km')
ylabel('gravity potential, J/kg')
title(' Solution of 1D Poisson equation')
