clear all
clc
%Datos
%La siguient línea está mal
disp('Modificado desde el mas halla...');
terms=50;
time=[0:0.01:5];
To=0.628;
syms t
pt1=15*t^0;
pt2=-15*t^0;
%Datos del SDF system
Tn=To/4;
k=3600;
dam=0.05;
%C�lculos
wo=2*pi/To;
ao=(1/To)*int(pt1,0,To/2)+(1/To)*int(pt2,To/2,To);

    aj=zeros(1,length(terms));
    bj=zeros(1,length(terms));
    ajf=zeros(length(terms),length(time));
    bjf=zeros(length(terms),length(time));
    p=zeros(1,length(time));
    u=zeros(1,length(time));
for i=1:length(time)
    for j=1:terms
        aj(j)=(2/To)*int(pt1*cos(j*wo*t),0,To/2)+(2/To)*int(pt2*cos(j*wo*t),To/2,To);
        bj(j)=(2/To)*int(pt1*sin(j*wo*t),0,To/2)+(2/To)*int(pt2*sin(j*wo*t),To/2,To);
        ajf(j,i)=aj(j)*cos(j*wo*time(i));
        bjf(j,i)=bj(j)*sin(j*wo*time(i));
        %Respuesta del sistema de un solo grado de libertad
        bethaj=j*wo/k;
        term(j,i)=(1/k)*(1/((1-bethaj^2)^2+(2*dam*bethaj)^2))*(((aj(j)*(2*dam*bethaj)+bj(j)*(1-bethaj^2))*sin(j*wo*time(i)))+(aj(j)*(1-bethaj^2)-bj(j)*2*dam*bethaj)*cos(j*wo*time(i)));  
    end
    p(i)=double(ao+sum(ajf(:,i))+sum(bjf(:,i)));    
    u(i)=double((ao/k)+sum(term(:,i)));
end
subplot(2,1,1);
plot(time,p,'b')
title('Fuerza aplicada')
xlabel('Tiempo')
ylabel('Fuerza')
hold on
subplot(2,1,2);
plot(time,u,'b')
title('Respuesta del sistema')
xlabel('Tiempo')
ylabel('Desplazamiento')
hold on

