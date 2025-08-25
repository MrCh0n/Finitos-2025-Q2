function [K] = crearK_CST(nodos, C)
%Crea la matriz de resistencia de un CST sin el espesor
%nodos es una lista (x, y) de los 3 nodos del elemento
arguments (Input)
   nodos(3,2) {mustBeNumeric}
   C(3,3) {mustBeNumeric}
end

x = nodos(:,1);
y = nodos(:,2);
one = ones(size(nodos,1),1);
A = [one x y];

Bx = [0 1 0]/A;
By = [0 0 1]/A;

B(1, [1 3 5]) = Bx;
B(2, [2 4 6]) = By;
B(3, :) = [By(1) Bx(1) By(2) Bx(2) By(3) Bx(3)];

Area = norm(cross([nodos(2,:) 0],[nodos(3,:) 0]))/2;

K = Area*B'*C*B;
end