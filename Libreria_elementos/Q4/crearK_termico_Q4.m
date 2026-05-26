function [K] = crearK_termico_Q4(nodos,k)
%Crea la matriz de resistencia de un Q4 sin espesor
%
%K = crearK_Q4(nodos, k)
%
%nodos es una lista (x, y) de los 4 nodos del elemento
%k es la conductividad del material
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
A = inv(A);


%% Gauss
puntos = [-sqrt(3/5) 0 sqrt(3/5)];
w = [5/9 8/9 5/9];

orden = size(puntos,2);

B = zeros(2,cant_puntos);
K = 0;
%los del for loop precalculado le saca la mitad de la corrida al crearK

% mat_Neta = [-0.4436    0.4436    0.0564   -0.0564;  
%             -0.2500    0.2500    0.2500   -0.2500;
%             -0.0564    0.0564    0.4436   -0.4436];
% 
% mat_Nzeta = [ -0.4436   -0.0564    0.0564    0.4436;
%                -0.2500   -0.2500    0.2500    0.2500;
%                -0.0564   -0.4436    0.4436    0.0564];

for i = 1:orden
    for j = 1:orden
        Neta = [0, 1, 0 puntos(j)]*A;
        Nzeta = [0, 0, 1, puntos(i)]*A;

        % Neta = mat_Neta(j,:);
        % Nzeta = mat_Nzeta(i,:);
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
    
        Bs = J\D;
            
        Bx = Bs(1,:);
        By = Bs(2,:);
        

        %crear la matrz B
        B(1,:) = Bx;
        B(2,:) = By;

        mult = abs(det(J))*w(i)*w(j);
        Kmin = B'*k*B;
        
        K = K + Kmin*mult;

    end% j
end% i

end