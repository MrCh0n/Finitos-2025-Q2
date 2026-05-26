function [M] = global_Q8(nodos)
%Crea la matriz de NT*N de un Q8 para el suavizado global o patch
%
%M = global_Q8(nodos)
%
% nodos es una lista (x, y) de los 4 nodos del elemento
%
%poner los nodos como:
% 4 - 7 - 3
% |       |
% 8       6
% |       |
% 1 - 5 - 2 

arguments (Input)
   nodos(8,2) {mustBeNumeric}
end
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
[w, puntos, n] = gauss([3,3]);

M=0;
for i = 1:n
        xi = puntos(i,1);
        eta = puntos(i,2);
    
        N = [1, xi, eta, xi^2, xi*eta,  eta^2,  xi^2*eta,  xi*eta^2]*A;

        % Neta = [0, 1, 0 puntos(i,2)]*A;
        % Nzeta = [0, 0, 1, puntos(i,1)]*A;
        % 
        % D = [Neta; Nzeta];
        % 
        % J = D*nodos;
  
        mult = w(i);
        Mmin = N'*N;
        
        M = M + Mmin*mult;
end% i

end