clc;
clear;

valores = struct;

Coord = [0 0;
    1 0;
    2 0;
    0 1;
    1 1;
    2 1]; %coordenadas x y de los nodos
%Coord = [0 0; 1 0; 0 1];


Elem = [1 2;
    2 3;
    1 4;
    2 4;
    2 5;
    2 6;
    3 6;
    4 5;
    5 6]; %referir a que nodos estan conectados los elementos
%Elem = [1 2; 2 3; 3 1];

valores.coord = Coord;
valores.elem = Elem;

valores
valores.coord(2, :)

valores = rmfield(valores, 'elem');

valores

syms x L
X = [1 x x^2 x^3];
A = [subs(X, x ,0)
     subs(diff(X,x), x, 0)
     subs(X, x, L)
     subs(diff(X,x), x, L)];

N = X/A;