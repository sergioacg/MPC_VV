function [sys,x0]=smodelo(t,x,u,flag,k10,k20,k30,e1,e2,e3,dhab,dhbc,dhad,ro,cp,kw,ar,ca0,T0,v,cae,cbe,Te)
%% Español
% En el Flag 0 es el paso donde indicamos para la S-Function que es lo que
% la función va a encontrar en el momendo de leer el modelo. Esa
% información será encontrada en un vector de 6 elementos, que llamaremos
% "sys"
%% Portugues
%O Flag 0 é o passo onde indicamos para a S-Function o que ela vai
%encontrar quando ler o modelo. Isso é feito com un vetor de 6 elementos
%chamado nesse caso como sys
if flag==0
   %% Español
   %Elemento 1: Número de estados continuos (Ecuaciones diferenciables)= 2
   %Elemento 2: Número de estados discretos: 0
   %Elemento 3: Número de saídas do modelo: 2
   %Elemento 4: Número de Entradas do modelo (a1) (a2): 2 entradas
   %Elemento 5: Parametro de control, colcar 0
   %Elemento 6: Tipologia do processo, (1 para processos continuos)
   
   %% Portugues
   %Elemento 1: quantidade de estados continuos (Equações diferenciaveis)= 3
   %Elemento 2: Numero de estados discretos: 0
   %Elemento 3: Numero de saídas do modelo: 3
   %Elemento 4: Numero de Entradas do modelo (FV) (Temp): 2 entradas
   %Elemento 5: Parametro de control, colcar 0
   %Elemento 6: Tipologia do processo, (1 para processos continuos)
   [sys]=[3,0,3,2,0,0];
   %Incluir condições iniciais
   x0=[cae cbe Te];
elseif flag==1
    %% Español Flag 1 llama el modelo
    %% Portugues Flag 1 chama o modelo
   [sys]=VanVusseModel(t,x,u,k10,k20,k30,e1,e2,e3,dhab,dhbc,dhad,ro,cp,kw,ar,ca0,T0,v);
elseif flag==3
    %% Español
    %Flag 3 indica la respuesta que se debe obtener, en este caso, un vetor con las 2
    %variables de estado
    
    %% Portugues
    %Flag 3 indica a resposta que deve obter, no caso, um vetor com as 3
    %variáveis de estado
   [sys]=x;
else
    %Como passo final com um vetor nulo é fechado o bucle
   [sys]=[];
end
end