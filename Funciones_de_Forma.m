clc
clear

syms x y x1 x2 x3 y1 y2 y3

%% Triangulo cte
nodo1 = [x1 y1];
nodo2 = [x2 y2];
nodo3 = [x3 y3];

X = [1; x; y];

A = [1 nodo1;
    1 nodo2;
    1 nodo3];

N = A\X

triangulo = [nodo1;nodo2;nodo3]

plot(polyshape(triangulo), 'Facecolor', 'g')

%% Cuadrado
%Q4
nodo1 = [0.5 0.5];
nodo2 = [-0.5 0.5];
nodo3 = [0.5 -0.5];
nodo4 = [-0.5 -0.5];

X = [1 x y x*y];

A = [subs(X, [x y], nodo1);
    1 nodo2 nodo2(1)*nodo2(2);
    1 nodo3 nodo3(1)*nodo3(2);
    1 nodo4 nodo4(1)*nodo4(2)];

N = X/A