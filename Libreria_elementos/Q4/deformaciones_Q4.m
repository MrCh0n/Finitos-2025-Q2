function [epsilon] = deformaciones_Q4(nodos,Uel)
%Devuelve un vector [exx1, eyy1, exy1, exx2..., exy4] de deformaciones de un Q4
%
%epsilon = deformaciones_Q4(nodos, Uel)
%
%nodos es una lista (x, y) de los 4 nodos del elemento
%U es un vector [u1, v1, u2..., v4] de los desplazamientos de los nodos
%
%poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 
%% creo el isoparametrico
cant_puntos = 4;

x1 = [-1; 1; 1; -1];
y1 = [-1; -1; 1; 1];
A = [ones(cant_puntos,1) x1 y1 x1.*y1];
A=inv(A);


%% Gauss
[w, puntos, n] = gauss([1,1]);

Bel = zeros(3,2*cant_puntos);
 
dir1 = 1:2:2*cant_puntos;
dir2 = 2:2:2*cant_puntos;
epsilon = zeros(3*cant_puntos,1);

for i = 1:n
    Neta = [0, 1, 0 puntos(i,2)]*A;
    Nzeta = [0, 0, 1, puntos(i,1)]*A;
    
    D = [Neta; Nzeta];

    J = D*nodos;

    Bs = J\D;
        
    Bx = Bs(1,:);
    By = Bs(2,:);

    %crear la matrz B
    Bel(1,dir1) = Bx;
    Bel(2,dir2) = By;
    Bel(3,dir1) = By;
    Bel(3,dir2) = Bx;
    
    e(i,:) = Bel*Uel;%epsilon del punto de gauss
    xy = e([1 2 3])';
    epsilon = epsilon + [xy; xy; xy; xy];%TODO deberia extrapolar a las esquinas pero es en (0,0) por ahora
end% i 
end