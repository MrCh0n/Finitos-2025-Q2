function [K] = crearK_H8(nodos,C)
%Crea la matriz de resistencia de un H8
%
%K = crearK_H8(nodos, C)
%
%nodos es una lista (x, y) de los 8 nodos del elemento
%C es el constitutivo de 6x6
%
%poner los nodos como:
%
%    8-------7
%  / |     / |
% 5 ------6  |
% |  |    |  |
% |  4----|--3 
% | /     | / 
% 1 - - - 2 
%          z ^   y
% con ejes   | ↗
%             --> X

arguments (Input)
   nodos(8,3) {mustBeNumeric}
   C(6,6) {mustBeNumeric}
end
%% creo el isoparametrico
cant_puntos = 8;

x1 = [-1; 1; 1; -1; -1; 1; 1; -1];
y1 = [-1; -1; 1; 1; -1; -1; 1; 1];
z1 = [-1; -1; -1; -1; 1; 1; 1; 1];
A = [ones(cant_puntos,1) x1 y1 z1 x1.*y1 x1.*z1 y1.*z1 x1.*z1.*y1];


%% Gauss
[w, puntos, n] = gauss([3,3]);

B = zeros(6,3*cant_puntos);
K = 0;
A = inv(A);

for i = 1:n
    %estan en los naturales todavia
    %pero para no confundirme
    Nx = [0, 1, 0, 0, puntos(i,2), puntos(i,3), 0, puntos(i,2)*puntos(i,3)]*A;
    Ny = [0, 0, 1, 0, puntos(i,1), 0, puntos(i,3), puntos(i,3)*puntos(i,1)]*A;
    Nz = [0, 0, 0, 1, 0, puntos(i,1), puntos(i,2), puntos(i,2)*puntos(i,1)]*A;

    D = [Nx; Ny; Nz];

    J = D*nodos;
    
    %paso las derivadas a x,y,z
    Bs = J\D;
        
    Bx = Bs(1,:);
    By = Bs(2,:);
    Bz = Bs(2,:);
    
    dirX = 1:3:3*cant_puntos;
    dirY = 2:3:3*cant_puntos;
    dirZ = 3:3:3*cant_puntos;

    %crear la matrz B
    B(1,dirX) = Bx;
    B(2,dirY) = By;
    B(3,dirZ) = Bz;

    B(4,dirY) = Bx;
    B(4,dirX) = By;

    B(5,dirZ) = Bx;
    B(5,dirX) = Bz;

    B(6,dirY) = Bz;
    B(6,dirZ) = By;

    mult = abs(det(J))*w(i);
    Kmin = B'*C*B;
    
    K = K + Kmin*mult;
end% i

end
