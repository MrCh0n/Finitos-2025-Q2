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

%A sin syms
x1 = nodos(:,1);
y1 = nodos(:,2);
one = ones(size(nodos,1),1);
A = [one x1 y1 x1.^2 x1.*y1 y1.^2 x1.^2.*y1 y1.^2.*x1];

A = inv(A);
%% GAUSS
tic
puntos = [-sqrt(3/5) 0 sqrt(3/5)];
w = [5/9 8/9 5/9];

x2 = (puntos+1)*(xmax-xmin)/2 + xmin;
y2 = (puntos+1)*(ymax-ymin)/2 + ymin;
integral = 0;
for i = 1:3
    for j = 1:3
        Bx = [0, 1, 0, 2*x2(i), y2(j), 0, 2*x2(i)*y2(j), y2(j)^2]*A;
        By = [0, 0, 1, 0, x2(i), 2*y2(j), x2(i)^2, 2*x2(i)*y2(j)]*A;
        
        dir1 = 1:2:16;
        dir2 = 2:2:16;
        B(1,dir1) = Bx;
        B(2,dir2) = By;
        B(3,dir1) = By;
        B(3,dir2) = Bx;

        integral = integral + B'*C*B*w(i)*w(j);
    end%for j
end%for i
toc

%% Syms
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
K = eval(K);
(integral-K)./K
end