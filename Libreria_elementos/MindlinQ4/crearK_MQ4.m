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

%Bending
puntos = [-sqrt(1/3) sqrt(1/3)];
w = [1 1];

orden = size(puntos,2);

for i = 1:orden
    for j = 1:orden
        Neta = [0, 1, 0 puntos(j)]/A;
        Nzeta = [0, 0, 1, puntos(i)]/A;
        
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
    
        mult = abs(det(J))*w(i)*w(j);
        Kmin_b = B_b'*D_b*B_b;
        
        K_b = K_b + Kmin_b*mult;

    end% j
end% i


%Shear
puntos = [0];
w = [2];
%puntos = [-sqrt(1/3) sqrt(1/3)];
%w = [1 1];

orden = size(puntos,2);
for i = 1:orden
    for j = 1:orden
      
        N = [1,puntos(i),puntos(j),puntos(i)*puntos(j)]/A;
        Neta = [0, 1, 0 puntos(j)]/A;
        Nzeta = [0, 0, 1, puntos(i)]/A;
        
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
    
        mult = abs(det(J))*w(i)*w(j);
        Kmin_s = B_s'*D_s*B_s;
        
        K_s = K_s + Kmin_s*mult;

    end% j
end% i

K = K_b + K_s;
end