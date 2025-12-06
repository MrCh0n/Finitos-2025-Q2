function [K] = crearK_termico_CST(nodos, k)
%Crea la matriz de resistencia local de un CST sin el espesor
%
%K = crearK_CST(nodos, C)
%
%nodos es una lista (x, y) de los 3 nodos del elemento
%k es la contuctividad termica
%
%poner los nodos como:
%    3
%   / \
%  /   \
% /     \
%1 - - - 2 

%% creo el isoparametrico
cant_puntos = 3;

x1 = [0; 1; 0];
y1 = [0; 0; 1];
A = [ones(cant_puntos,1) x1 y1];


%% Gauss
puntos =     [1/6 2/3
    1/6    1/6
    2/3    1/6];

w = [1 1 1]/6;

orden = size(puntos,1);

B = zeros(2,cant_puntos);
K = 0;

for i = 1:orden
    Neta = [0, 1, 0]/A;
    Nzeta = [0, 0, 1]/A;
    
    D = [Neta; Nzeta];

    J = D*nodos;

    Bs = J\D;
        
    Bx = Bs(1,:);
    By = Bs(2,:);
    
    %crear la matrz B
    B(1,:) = Bx;
    B(2,:) = By;
    
    mult = abs(det(J))*w(i);
    Kmin = B'*k*B;
  
    K = K + Kmin*mult;
end

end