function dydt = VanVusseModel(t,x,u,k10,k20,k30,E1,E2,E3,deltaAB,deltaBC,deltaAD,ro,cp,Kw,Ar,Ca0,T0,V)
%% Función que representa la dinámica del reactor de Van de Vusse
% x = Estados del sistema
% u = Estructura con las dos entradas del sistema

%% Datos   
% Entradas
fov = u(1);          %h^-1  (Flujo)
Tk  = u(2);      %ºC    (Temperatura de la camisa)

% 
%  Notación para las Variables de Estado
%
   ca   = x(1);     %Concentración A
   cb   = x(2);     %Concentración B
   T    = x(3);     %Temperatura
   
%Constantes cinéticas
k1 = k10*exp(E1./(T+273.15));
k2 = k20*exp(E2./(T+273.15));
k3 = k30*exp(E3./(T+273.15));

%
%  Ecuaciones Diferenciales del Sistema
%
dydt = zeros(3,1);
dydt(1)  = fov*(Ca0 - ca) -k1*ca - k3*ca*ca;
dydt(2)  = -fov*cb + k1*ca - k2*cb;
dydt(3)  = (1/(ro*cp))*(k1*ca*deltaAB + k2*cb*deltaBC + k3*ca^2*deltaAD)...
              + fov*(T0-T) + (Kw*Ar/(ro*cp*V))*(Tk-T);

