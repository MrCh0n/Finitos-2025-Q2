function [K] = crearK_Q4(nodos,C)
%Crea la matriz de resistencia de un Q4 sin espesor
%
%K = crearK_Q4(nodos, C)
%
%nodos es una lista (x, y) de los 4 nodos del elemento
%C es el constitutivo de 3x3
%
%poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 

arguments (Input)
   nodos(4,2) {mustBeNumeric}
   C(3,3) {mustBeNumeric}
end
%% creo el isoparametrico
cant_puntos = 4;

x1 = [-1; 1; 1; -1];
y1 = [-1; -1; 1; 1];
A = [ones(cant_puntos,1) x1 y1 x1.*y1];


%% Gauss
puntos = [-sqrt(3/5) 0 sqrt(3/5)];
w = [5/9 8/9 5/9];

orden = size(puntos,2);

B = zeros(3,2*cant_puntos);
K = 0;

for i = 1:orden
    for j = 1:orden
        Neta = [0, 1, 0 puntos(j)]/A;
        Nzeta = [0, 0, 1, puntos(i)]/A;
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
    
        Bs = J\D;
            
        Bx = Bs(1,:);
        By = Bs(2,:);
        
        dir1 = 1:2:2*cant_puntos;
        dir2 = 2:2:2*cant_puntos;
        %crear la matrz B
        B(1,dir1) = Bx;
        B(2,dir2) = By;
        B(3,dir1) = By;
        B(3,dir2) = Bx;
    
        mult = abs(det(J))*w(i)*w(j);
        Kmin = B'*C*B;
        
        K = K + Kmin*mult;

    end% j
end% i

end