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
alfa = 1e-4;
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

function [K] = crearK_Q4(nodos,C)
%Crea la matriz de resistencia de un Q4 sin espesor
%
%K = crearK_Q4(nodos, C)
%
%nodos es una lista (x, y) de los 4 nodos del elemento
%C es el constitutivo de 3x3
%
%poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 

arguments (Input)
   nodos(4,2) {mustBeNumeric}
   C(3,3) {mustBeNumeric}
end
%% creo el isoparametrico
cant_puntos = 4;

x1 = [-1; 1; 1; -1];
y1 = [-1; -1; 1; 1];
A = [ones(cant_puntos,1) x1 y1 x1.*y1];


%% Gauss
[w, puntos, n] = gauss([2,2]);

B = zeros(3,2*cant_puntos);
K = 0;

for i = 1:n
        Neta = [0, 1, 0 puntos(i,2)]/A;
        Nzeta = [0, 0, 1, puntos(i,1)]/A;
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
    
        Bs = J\D;
            
        Bx = Bs(1,:);
        By = Bs(2,:);
        
        dir1 = 1:2:2*cant_puntos;
        dir2 = 2:2:2*cant_puntos;
        %crear la matrz B
        B(1,dir1) = Bx;
        B(2,dir2) = By;
        B(3,dir1) = By;
        B(3,dir2) = Bx;
    
        mult = abs(det(J))*w(i);
        Kmin = B'*C*B;
        
        K = K + Kmin*mult;
end% i

end

function [K] = crearK_MQ4(nodos,E,v,t)
%Crea la matriz de resistencia de una placa Mindlin Q4
%
%K = crearK_Q4(nodos, C)
%
%nodos es una lista (x, y) de los 4 nodos del elemento
%C es el constitutivo de 3x3
%
%poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 

arguments (Input)
   nodos(4,2) {mustBeNumeric}
   E {mustBeNumeric}
   v {mustBeNumeric}
   t {mustBeNumeric}
end
%% matriz constitutiva
G = E/(2*(1+v));
D = E*t^3/(12*(1-v^2));
k = 5/6;

D_b = [D v*D 0;
      v*D D 0;
      0 0 (1-v)/2*D];
D_s = [k*G*t 0;
      0 k*G*t];

%% creo el isoparametrico
cant_puntos = 4;

x1 = [-1; 1; 1; -1];
y1 = [-1; -1; 1; 1];
A = [ones(cant_puntos,1) x1 y1 x1.*y1];


%% Gauss
%puntos = [-sqrt(3/5) 0 sqrt(3/5)];
%w = [5/9 8/9 5/9];

B_b = zeros(3,3*cant_puntos);
B_s = zeros(2,3*cant_puntos);

dir1 = 1:3:3*cant_puntos;
dir2 = 2:3:3*cant_puntos;
dir3 = 3:3:3*cant_puntos;
K_b = 0;
K_s = 0;

%Bending (full)
[w, puntos, n] = gauss([2,2]);

for i = 1:n
        Neta = [0, 1, 0 puntos(i,2)]/A;
        Nzeta = [0, 0, 1, puntos(i,1)]/A;
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
    
        Bs = J\D;
            
        Bx = Bs(1,:);
        By = Bs(2,:);
        
        %dir1 = 1:2:3*cant_puntos;
        %dir2 = 2:2:3*cant_puntos;
        %dir3 = 3:2:3*cant_puntos; definidas afuera del loop
        %crear la matrz B_b
        B_b(1,dir2) = Bx;
        B_b(2,dir3) = By;
        B_b(3,dir2) = By;
        B_b(3,dir3) = Bx;
    
        mult = abs(det(J))*w(i);
        Kmin_b = B_b'*D_b*B_b;
        
        K_b = K_b + Kmin_b*mult;
end% i


%Shear 1x1
[w, puntos, n] = gauss([1,1]);

for i = 1:n
        N = [1,puntos(i,1),puntos(i,2),puntos(i,1)*puntos(i,2)]/A;
        Neta = [0, 1, 0 puntos(i,2)]/A;
        Nzeta = [0, 0, 1, puntos(i,1)]/A;
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
    
        Bs = J\D;
            
        Bx = Bs(1,:);
        By = Bs(2,:);
        
        %dir1 = 1:2:3*cant_puntos;
        %dir2 = 2:2:3*cant_puntos;
        %dir3 = 3:2:3*cant_puntos; definidas afuera del loop
        %crear la matrz B_b
        B_s(1,dir1) = -Bx;
        B_s(2,dir1) = -By;
        B_s(1,dir2) = N;
        B_s(2,dir3) = N;
    
        mult = abs(det(J))*w(i);
        Kmin_s = B_s'*D_s*B_s;
        
        K_s = K_s + Kmin_s*mult;
end% i

K = K_b + K_s;
end