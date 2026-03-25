function [K] = crearK_shellMQ4(nodos,E,v,t)
%Crea la matriz de resistencia de una cascara Mindlin Q4
%
%K = crearK_Q4(nodos, C)
%
%nodos es una lista (x, y, z) de los 4 nodos del elemento
%E es el modulo de elasticidad
%v el coef de poisson
%t el espesor
%
%poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 

arguments (Input)
   nodos(4,3) {mustBeNumeric}
   E {mustBeNumeric}
   v {mustBeNumeric}
   t {mustBeNumeric}
end

addpath(pwd+"/../MindlinQ4")
addpath(pwd+"/../Q4")

%plane stress
C = E*t/(1-v^2)*[1 v 0;
               v 1 0;
               0 0 (1-v)/2];

cant = 4; %cant de nodos x elem

%% Matriz de rotacion
v1 = (nodos(2,:)-nodos(1,:))/norm(nodos(2,:)-nodos(1,:));

v3 = cross(v1,(nodos(4,:)-nodos(1,:)));

v3 = v3 / norm(v3); % Normalize the cross product vector

v2 = cross(v3,v1);

R = [v1;v2;v3]; %v1 es vector fila

%% Rotar nodos
nodosp = (R*nodos')';

%desprecio zp
KQ4 = crearK_Q4(nodosp(:,1:2),C);
Kmindlin = crearK_MQ4(nodosp(:,1:2),E,v,t);

Kp = zeros(5*cant);

plano_xy = [1:2 6:7 11:12 16:17];
z_giros = [3:5 8:10 13:15 18:20];

Kp(plano_xy, plano_xy) = KQ4;
Kp(z_giros, z_giros) = Kmindlin;

%spy(Kp)

%% Rotar matriz
L = zeros(5,6);
L(1:3,1:3) = R;
L(4:5,4:6) = [-v2;v1];

T = zeros(5*cant,6*cant);
for i=1:cant
    T((i-1)*5+1:i*5, (i-1)*6+1:i*6) = L;
end
K = T'*Kp*T;


end