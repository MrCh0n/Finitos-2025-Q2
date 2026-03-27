clc
clear

%% Datos problema
E = 200e9;
t = 1;
v = 0.3;
%Plain stress
C = t*E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];

a = 1;

bordes = [0 0;
          a 0;
          a a;
          0 a];

divx = 10;
divy = 10;

%% creo el mesh y calculo la K
mesh = Q4(bordes,divx,divy,C);

mesh.armar_K;

%% hacer la carga
carga = [0 0 0 0;
         0 0  0 0;
         0 5860 0 5860;
         0 0 0 0]*0;

volumen = [zeros(4,1) -ones(4,1)];
% i puede ser el elemento justo que se quiere cargar, puse todos porque
% quiero hacer una fuerza volumetrica solo
for i = 1:mesh.counts.nelem
   mesh.armar_R(i,carga,volumen,true);
end

%% Condiciones de borde
% talvez se puede hacer mas facil, es empotrados en la parte de abajo
borde_inferior = mesh.nodes(:,2)==0;
for i = 1:mesh.counts.nnod
   if borde_inferior(i)
    mesh.cond_borde(i,3);
   end
end

%% Calculo y plot
mesh.calc_U;

%Dibuja ~5% de escala del modelo
mesh.dibujar;
