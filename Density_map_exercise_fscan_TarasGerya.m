clear all;close all;
%% extraction des données liés au manteau
fdata=fopen('m895_ro','rt');%lecture du fichier txt 'rt' pour read txt
a=fscanf(fdata,'%s',1);%lecture du 1er mot '%s' pour string et exclusion automatique passage au deuxième mot
%lecture des données floting point (nombre à virgule)'%f'
tnum=fscanf(fdata,'%f',1);% nombre de point dans la grille y (350x350)pour axe T
pnum=fscanf(fdata,'%f',1);% nombre de point dans la grille x (350x350)pour axe P
Tmin=fscanf(fdata,'%f',1);% min de température
Pmin=fscanf(fdata,'%f',1);% minimum de pression
Tstep=fscanf(fdata,'%f',1);%pas de Temperature
Pstep=fscanf(fdata,'%f',1);%pas de pression

b=fscanf(fdata,'%s',2);%exclusion de ce string 

rhomantle=zeros(pnum,tnum); %definition de la dimension de la matrice de densité 
for i=1:1:pnum
    for j=1:1:tnum
    rhomantle(i,j)=fscanf(fdata,'%f',1); %definition des valeurs de densité
    end
end

fclose(fdata); %fermeture du fichier avant ouverture d'un autre

%Même chose mais pour la croûte
fdata=fopen('morn_ro','rt');
a2=fscanf(fdata,'%s',1);
tnum2=fscanf(fdata,'%f',1);
pnum2=fscanf(fdata,'%f',1);
Tmin2=fscanf(fdata,'%f',1);
Pmin2=fscanf(fdata,'%f',1);
Tstep2=fscanf(fdata,'%f',1);
Pstep2=fscanf(fdata,'%f',1);

b2=fscanf(fdata,'%s',2);
rhocrust=zeros(pnum2,tnum2);
for i=1:1:pnum2
    for j=1:1:tnum2
    rhocrust(i,j)=fscanf(fdata,'%f',1);
    end
end
fclose(fdata);

%calcule de la pression et de la température 
figure(1)
for i=1:1:pnum
    P(i)=Pmin*1e-4+(i-1)*Pstep*1e-4;%Pression(i) GPa
end
for j=1:1:tnum
    T(j)=Tmin-273+(j-1)*Tstep; %Temperature(i) en °C
end

subplot(2,2,1);pcolor(T,P,rhomantle);title('Pyrolite density, kg/m^3');
xlabel('T,C');ylabel('P,GPa');shading interp;colorbar;colormap(jet);set(gca,'Ydir','reverse');
subplot(2,2,2);pcolor(T,P,rhocrust);title('MORB density, kg/m^3');
xlabel('T,C');ylabel('P,GPa');shading interp;colorbar;colormap(jet);set(gca,'Ydir','reverse');
subplot(2,2,3);pcolor(T,P,rhomantle-rhocrust);title('Pyrolite-MORB density difference, kg/m^3');
xlabel('T,C');ylabel('P,GPa');shading interp;colorbar;colormap(jet);set(gca,'Ydir','reverse');


