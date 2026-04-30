function [M] = global_Q4(nodos)
%Crea la matriz de NT*N de un Q4 para el suavizado global o patch
%
%M = global_Q4(nodos)
%
% nodos es una lista (x, y) de los 4 nodos del elemento
%
%poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 

arguments (Input)
   nodos(4,2) {mustBeNumeric}
end
%% creo el isoparametrico
cant_puntos = 4;

x1 = [-1; 1; 1; -1];
y1 = [-1; -1; 1; 1];
A = [ones(cant_puntos,1) x1 y1 x1.*y1];
A=inv(A);
%% Gauss
[w, puntos, n] = gauss([2,2]);

M=0;
for i = 1:n
        N = [1, puntos(i,1), puntos(i,2), puntos(i,1)*puntos(i,2)]*A;

        Neta = [0, 1, 0 puntos(i,2)]*A;
        Nzeta = [0, 0, 1, puntos(i,1)]*A;

        D = [Neta; Nzeta];

        J = D*nodos;
  
        mult = abs(det(J))*w(i);
        Mmin = N'*N;
        
        M = M + Mmin*mult;
end% i

end