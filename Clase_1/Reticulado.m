clear
clc

valores = struct;

%coordenadas x y de los nodos
valores.coord = [0 0;
    1 0;
    2 0;
    0 1;
    1 1;
    2 1]; 

%referir a que nodos estan conectados los elementos
valores.elem = [1 2;
    2 3;
    1 4;
    2 4;
    2 5;
    2 6;
    3 6;
    4 5;
    5 6]; 

%Condiciones de borde
% nodos agarrados (cada nodo tiene 2 elemento)
valores.BC = [1 7 8];

%Fuerzas externas
valores.Fnodos = [4 -2000;
                 6 -2000];

barras(valores);