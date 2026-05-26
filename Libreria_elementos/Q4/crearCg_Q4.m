function [Cg] = crearCg_Q4(nodos,t, alpha)
%Crea la matriz de acople de un Q4
%
%Cg = crearCg_Q4(nodos, t, aplha)
%
%nodos es una lista (x, y) de los 4 nodos del elemento
%t es el espesor del elemento
%alpha es el coeficiente de Biot
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
[w, puntos, n] = gauss([2,2]);

B = zeros(3,2*cant_puntos);
Cg = 0;
dir1 = 1:2:2*cant_puntos;
dir2 = 2:2:2*cant_puntos;

m = [1; 1; 0];

for i = 1:n
        xi = puntos(i,1);
        eta = puntos(i,2);

        N = [1, xi, eta, xi*eta]*A;

        Neta = [0, 1, 0, eta]*A;
        Nzeta = [0, 0, 1, xi]*A;
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
    
        Bs = J\D;
            
        Bx = Bs(1,:);
        By = Bs(2,:);

        %crear la matrz B
        B(1,dir1) = Bx;
        B(2,dir2) = By;
        B(3,dir1) = By;
        B(3,dir2) = Bx;
    
        mult = abs(det(J))*w(i)*t;
        Kmin = B'*alpha*m*N;
        
        Cg = Cg + Kmin*mult;
end% i

end