function [M] = masa_Q4(nodos,p,lumped)
%Crea la matriz de masa de un Q4 sin espesor
%
%K = masa_Q4(nodos, p, lumped)
%
% nodos es una lista (x, y) de los 4 nodos del elemento
% p es la densidad del elemento
% lumped devuelve la matriz completa si es 0 y la lumpeada si es otro
% numero
%
%poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 

arguments (Input)
   nodos(4,2) {mustBeNumeric}
   p(1,1) {mustBeNumeric}
   lumped(1,1) {mustBeNumeric}
end
%% creo el isoparametrico
cant_puntos = 4;

x1 = [-1; 1; 1; -1];
y1 = [-1; -1; 1; 1];
A = [ones(cant_puntos,1) x1 y1 x1.*y1];
A=inv(A);
%% lumpeado
if lumped
% Gauss
[w, puntos, n] = gauss([2,2]);

Area = 0;
for i = 1:n
        Neta = [0, 1, 0 puntos(i,2)]*A;
        Nzeta = [0, 0, 1, puntos(i,1)]*A;

        D = [Neta; Nzeta];

        J = D*nodos;
  
        mult = abs(det(J))*w(i);
        
        Area = Area + mult;
end% i
    m = Area*p/4;
    M = m*eye(2*cant_puntos);
    return
end


%% Gauss
[w, puntos, n] = gauss([2,2]);

M=0;
dir1 = 1:2:8;
dir2 = 2:2:8;
for i = 1:n
        N = [1, puntos(i,1), puntos(i,2), puntos(i,1)*puntos(i,2)]*A;

        Neta = [0, 1, 0 puntos(i,2)]*A;
        Nzeta = [0, 0, 1, puntos(i,1)]*A;

        D = [Neta; Nzeta];

        J = D*nodos;
    
        Ns(dir1,1) = N;
        Ns(dir2,2) = N;
  
        mult = abs(det(J))*w(i);
        Mmin = Ns*p*Ns';
        
        M = M + Mmin*mult;
end% i

end