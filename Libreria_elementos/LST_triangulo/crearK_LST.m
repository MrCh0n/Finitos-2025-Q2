function [K] = crearK_LST(nodos,C)
%Crea la matriz de resistencia local de un LST sin el espesor
%
%K = crearK_LST(nodos, C)
%
%nodos es una lista (x, y) de los 6 nodos del elemento
%C es el constitutivo de 3x3
%
%poner los nodos como:
%    3
%   / \
%  6   5
% /     \
%1 - 4 - 2 

arguments (Input)
   nodos(6,2) {mustBeNumeric}
   C(3,3) {mustBeNumeric}
end

%% creo el isoparametrico
cant_puntos = 6;

x1 = [0; 1; 0; 0.5; 0; 0.5];
y1 = [0; 0; 1; 0; 0.5; 0.5];
A = [ones(cant_puntos,1) x1 y1 x1.*y1 x1.^2 y1.^2];

%N = [1 x y x.*y x.^2 y.^2]/A;
%% Gauss
puntos =     [1/6 2/3
    1/6    1/6
    2/3    1/6];

w = [1 1 1]/6;

orden = size(puntos,1);

B = zeros(3,cant_puntos);
K = 0;

for i = 1:orden
    Neta = [0, 1, 0, puntos(i,2), 2*puntos(i,1), 0]/A;
    Nzeta = [0, 0, 1, puntos(i,1), 0, 2*puntos(i,2)]/A;
    
    J = [Neta; Nzeta]*nodos;%el jacobiano para el punto elegido
    
    Bx = Neta*J1(1,1) + Nzeta*J1(1,2);
    By = Neta*J1(2,1) + Nzeta*J1(2,2);
  
    dir1 = 1:2:2*cant_puntos;
    dir2 = 2:2:2*cant_puntos;
    %crear la matrz B
    B(1,dir1) = Bx;
    B(2,dir2) = By;
    B(3,dir1) = By;
    B(3,dir2) = Bx;
    
    mult = det(J)*w(i);
    Kmin = B'*C*B;

    K = K + Kmin*mult;
end

end