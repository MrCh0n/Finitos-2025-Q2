function [M] = masa_CST(nodos,p,lumped)
%Crea la matriz de masa de un CST sin espesor
%
%K = masa_LST(nodos, p, lumped)
%
% nodos es una lista (x, y) de los 3 nodos del elemento
% p es la densidad del elemento
% lumped devuelve la matriz completa si es 0 y la lumpeada si es otro
% numero
%
%poner los nodos como:
%    3
%   / \
%  /   \
% /     \
%1 - - - 2 

arguments (Input)
   nodos(3,2) {mustBeNumeric}
   p(1,1) {mustBeNumeric}
   lumped(1,1) {mustBeNumeric}
end
%% creo el isoparametrico
cant_puntos = 3;

x1 = [0; 1; 0];
y1 = [0; 0; 1];
A = [ones(cant_puntos,1) x1 y1];
A = inv(A);
%% lumpeado
if lumped
% Gauss
puntos =     [1/6 2/3
    1/6    1/6
    2/3    1/6];

w = [1 1 1]/6;
n=3;

Area = 0;
for i = 1:n
        Neta = [0, 1, 0]*A;
        Nzeta = [0, 0, 1]*A;

        D = [Neta; Nzeta];

        J = D*nodos;
  
        mult = abs(det(J))*w(i);
        
        Area = Area + mult;
end% i
    m = Area*p/3*eye(2*cant_puntos);
    M = m;
    return
end


%% Gauss
puntos =     [1/6 2/3
    1/6    1/6
    2/3    1/6];

w = [1 1 1]/6;

n=3;

M=0;
dir1 = 1:2:6;
dir2 = 2:2:6;
for i = 1:n
        N = [1, puntos(i,1), puntos(i,2)]*A;
     
        Neta = [0, 1, 0]*A;
        Nzeta = [0, 0, 1]*A;

        D = [Neta; Nzeta];

        J = D*nodos;
    
        Ns(dir1,1) = N;
        Ns(dir2,2) = N;
  
        mult = abs(det(J))*w(i);
        Mmin = Ns*p*Ns';
        
        M = M + Mmin*mult;
end% i

end