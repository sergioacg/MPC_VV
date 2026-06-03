%% Sergio Castaño
% Universidad de Pamplona - 2026
% http://controlautomaticoeducacion.com/
% Reactor de Van de Vusse NO isotermico proceso MIMO
%%
clc
close all
clear all

%Entradas
u(1)=72;            %Entrada do sistema 2 (Flujo) 72
u(2)=128.95;        %Entrada do sistema 1 (temp) 128.95

%Condiciones Iniciales
CA0d  = 5.1;          % mol/l
CB0d  = 0;            % mol/l
T0d   = 130;          % C
x0 = [CA0d CB0d T0d];

%Periodo de Muestreo
Ts=0.1;
nm=100;
t = 0:Ts:(nm-1)*Ts; %Tiempo de Simulación
tsim = 0:Ts:0.1;

%% Linealizacion del sistema por Jacobiana
%DADOS        
k10 = 1.287e12; %h^-1
k20 = 1.287e12; %h^-1
k30 = 9.043e9;  %L/molA.h
E1 = -9758.3; %K
E2 = -9758.3; %K
E3 = -8560.0; %K
deltaAB = -4.20; %kJ/molA
deltaBC = 11.00; %kJ/molB
deltaAD = 41.85; %kJ/molA
ro = 0.9342; %Kg/L
cp = 3.01; %Kj/Kg.K
Kw = 4032; %Kj/h.K.m^2
Ar = 0.215; %m^2
V = 10; %l

%Valores iniciales
T0 = x0(3); %ºC
Ca0 = x0(1); %molA/L

%Entradas
FV = u(1); %h^-1
Tk = u(2); %ºC

par = [k10 k20 k30 E1 E2 E3 deltaAB deltaBC deltaAD ro cp Kw Ar V Tk T0 Ca0 FV];


% Encuentra el Estado Estacionario para las entradas dadas
X = fsolve(@(x)VanVusseModel(t,x,u,k10,k20,k30,E1,E2,E3,deltaAB,deltaBC,deltaAD,ro,cp,Kw,Ar,Ca0,T0,V),x0);

%Establesco los Estados en el Estado Estacionario
Cae = X(1)
Cbe = X(2)
Te = X(3)

%Cálculo das constantes cinéticas no E.E.
K1e = k10*exp(E1/(Te+273.15));
K2e = k20*exp(E2/(Te+273.15));
K3e = k30*exp(E3/(Te+273.15));


%Cálculo da matriz A
a11 = -FV-K1e-2*K3e*Cae;
a12 = 0;
a13 = K1e*Cae*E1./(Te+273.15)^2 + K3e*(Cae^2)*E3./(Te+273.15)^2;

a21 = K1e;
a22 = -FV-K2e;
a23 = -K1e*Cae*E1./(Te+273.15)^2 + K2e*Cbe*E2./(Te+273.15)^2;

a31 = (1/(ro*cp))*(K1e*deltaAB+2*K3e*Cae*deltaAD);
a32 = (1/(ro*cp))*(K2e*deltaBC);
a33 = (1/(ro*cp))*(-K1e*Cae*deltaAB*E1./(Te+273.15)^2 - K2e*Cbe*deltaBC*E2./(Te+273.15)^2 - K3e*(Cae^2)*deltaAD*E3./(Te+273.15)^2)-FV-(Kw*Ar)/(ro*cp*V);

A = [a11 a12 a13     
     a21 a22 a23                                             
     a31 a32 a33];      

%Cálculo da matriz B
B = [Ca0-Cae  0
     -Cbe     0
     T0-Te   (Kw*Ar)/(ro*cp*V)];

%Cálculo da matriz C
C = [0 0 0
     0 1 0
     0 0 1];

%Cálculo da matriz D
D = [0 0
     0 0
     0 0];

% %Determino la Funcion de Transferencia para la entrada 1 Cb
% display('Representación en Función de Transferencia para Cb')
% [num1,den1]=ss2tf(A,B,C,D,1);
% ft1=tf(num1,den1)

H = ss(A,B,C,D);
[z,p,k] = zpkdata(H);
K = zpk(z,p,k)

[n1,d1] = ss2tf(A,B,C,D,1); %FT de y1/u1 e y2/u1
[n2,d2] = ss2tf(A,B,C,D,2); %FT de y1/u2 e y2/u2

G11=tf(n1(2,:),d1);
G12=tf(n2(2,:),d2);
G21=tf(n1(3,:),d1);
G22=tf(n2(3,:),d2);

Ps=[G11 G12;G21 G22];

%Comprobar Modelo Lineal VS Modelo NO Lineal
du=5;
u(1)=u(1)+du;
%Modo antiguio
% [t,x] = ode45(@(t,x)VanVusseModel(t,x,u,k10,k20,k30,E1,E2,E3,deltaAB,deltaBC,deltaAD,ro,cp,Kw,Ar,Ca0,T0,V), tsim , X, u);
% Modo Nuevo
[t,x] = ode45(@(t,x)VanVusseModel(t,x,u,k10,k20,k30,E1,E2,E3,deltaAB,deltaBC,deltaAD,ro,cp,Kw,Ar,Ca0,T0,V), tsim, X);

tsm=0:0.001:0.1;
in=zeros(length(tsm),2);
in(:,1)=du*ones(length(tsm),1);
in(:,2)=0*ones(length(tsm),1);
[ylin]=lsim(H,in,tsm);

%Grafica
figure
plot(t,x(:,2),'-b',tsm,ylin(:,2)+X(2),'--r','LineWidth',2),ylabel('C_b');
title('Comparación modelo Lineal vs No Lineal (C_b)');
xlabel('tiempo (min)');
ylabel('C_b mol/litro');
legend('No Lineal','Lineal','Location','SouthEast')
figure
plot(t,x(:,3),'-b',tsm,ylin(:,3)+X(3),'--r','LineWidth',2),ylabel('Q');
title('Comparación modelo Lineal vs No Lineal (T)');
xlabel('tiempo (min)');
ylabel('T (K)');
legend('No Lineal','Lineal','Location','SouthEast')