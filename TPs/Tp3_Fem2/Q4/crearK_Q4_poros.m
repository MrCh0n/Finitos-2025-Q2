function [K] = crearK_Q4_poros(nodos, k, mu, t)
%Crea la matriz de difusion de un Q4
%
%K = crearK_Q4_poros(nodos, k, mu)
%
%nodos es una lista (x, y) de los 4 nodos del elemento
%k es la permiabilidad
%mu es la viscosidad
%t es el espesor
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

K = 0;

for i = 1:n
        xi = puntos(i,1);
        eta = puntos(i,2);

        Neta = [0, 1, 0, eta]*A;
        Nzeta = [0, 0, 1, xi]*A;
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
    
        Bs = J\D;
            
        mult = abs(det(J))*w(i)*t;
        Kmin = k/mu*(Bs'*Bs);
        
        K = K + Kmin'*mult;
end% i

end