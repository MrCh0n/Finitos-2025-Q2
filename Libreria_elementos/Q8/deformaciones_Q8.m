function [epsilon] = deformaciones_Q8(nodos,Uel)
%Devuelve un vector [exx1, eyy1, exy1, exx2..., exy8] de deformaciones de
%un Q8
%
%epsilon = deformaciones_Q8(nodos, Uel)
%
%nodos es una lista (x, y) de los 8 nodos del elemento
%U es un vector [u1, v1, u2..., v8] de los desplazamientos de los nodos
%
%poner los nodos como:
% 4 - 7 - 3
% |       |
% 8       6
% |       |
% 1 - 5 - 2 
%% creo el isoparametrico
cant_puntos = 8;

puntos_Q8 = [-1    -1;
              1    -1;
              1     1;
             -1     1;
              0    -1;
              1     0;
              0     1;
             -1     0];

x1 = puntos_Q8(:,1);%puntos del cuadrado
y1 = puntos_Q8(:,2);
A = [ones(cant_puntos,1) x1 y1 x1.^2 x1.*y1 y1.^2 x1.^2.*y1 y1.^2.*x1];
A = inv(A);

%% Gauss
[w, puntos, n] = gauss([2,2]);

Bel = zeros(3,2*cant_puntos);
 
dir1 = 1:2:2*cant_puntos;
dir2 = 2:2:2*cant_puntos;
epsilon = zeros(3*cant_puntos,1);

for i = 1:n
    xi = puntos(i,1);
    eta = puntos(i,2);
    
    Neta = [0, 1, 0, 2*xi, eta, 0, 2*xi*eta, eta^2]*A;%derivada de N en eta en los puntos de Gauss
    Nzeta = [0, 0, 1, 0, xi, 2*eta, xi^2, 2*xi*eta]*A;%derivada de N en zeta
    
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

    mult = abs(det(J))*w(i);
    
    e(i,:) = Bel*Uel*mult;%epsilon del punto de gauss
    xy = e([1 2 3])';
    epsilon = epsilon + [xy; xy; xy; xy; xy; xy; xy; xy];%TODO deberia extrapolar a las esquinas pero es en (0,0) por ahora
end% i 
end