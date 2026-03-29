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
%area
Ae = area(nodosp(:,1:2));
alfa = 0.5;
G = 0.5*E/(1 + v);
Ktitaz = alfa*Ae*G*t*eye(4);

Kp = zeros(6*cant);

plano_xy = [1:2 7:8 13:14 19:20];
z_giros = [3:5 9:11 15:17 21:23];
tita_z = 6:6:24;

Kp(plano_xy, plano_xy) = KQ4;
Kp(z_giros, z_giros) = Kmindlin;
Kp(tita_z,tita_z) = Ktitaz;


%spy(Kp)

%% Rotar matriz
L = zeros(5,6);
L(1:3,1:3) = R;
L(4:6,4:6) = [-v2;v1;v3];

T = zeros(6*cant);
for i=1:cant
    T((i-1)*6+1:i*6, (i-1)*6+1:i*6) = L;
end
K = T'*Kp*T;


end

function [Ae] = area(nodos)
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(4,1) x1 y1 x1.*y1];
    [w, puntos, n] = gauss([2,2]);
    Ae = 0; % area del elemento
    for i = 1:n
        Neta = [0, 1, 0 puntos(i,2)]/A;
        Nzeta = [0, 0, 1, puntos(i,1)]/A;
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
        mult = abs(det(J))*w(i);
        Ae = Ae + mult;
    end
end