function [K] = crearK_Q8(nodos, C)
%Crea la matriz de resistencia de un Q8
%nodos es una lista (x, y) de los 8 nodos del elemento
arguments (Input)
   nodos(8,2) {mustBeNumeric}
   C(3,3) {mustBeNumeric}
end
syms x y

xmin = min(nodos(:,1));
xmax = max(nodos(:,1));
ymin = min(nodos(:,2));
ymax = max(nodos(:,2));

X = [1 x y x^2 x*y y^2 x^2*y y^2*x];
Xx = diff(X,x);
Xy = diff(X,y);

A = zeros(8);
for i = 1:8
    A(i,:) = subs(X, [x y], nodos(i,:));
end
A = inv(A);

Bx = Xx*A;
By = Xy*A;

B=[];
for i = 1:8
    dx = Bx(i);
    dy = By(i);
    Bfila = [dx 0; 0 dy; dy dx];
    B = [B Bfila];
end

integrando = B'*C*B;

tic
K = int(int(integrando, x, xmin, xmax), y, ymin, ymax);
toc

end