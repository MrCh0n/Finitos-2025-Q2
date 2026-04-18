function [M] = masa_Q8(nodos,p,lumped)
%Crea la matriz de masa de un Q4 sin espesor
%
%K = masa_Q4(nodos, p, lumped)
%
% nodos es una lista (x, y) de los 8 nodos del elemento
% p es la densidad del elemento
% lumped devuelve la matriz completa si es 0 y la lumpeada si es otro
% numero
%
%poner los nodos como:
% 4 - 7 - 3
% |       |
% 8       6
% |       |
% 1 - 5 - 2 

arguments (Input)
   nodos(8,2) {mustBeNumeric}
   p(1,1) {mustBeNumeric}
   lumped(1,1) {mustBeNumeric}
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
%% lumpeado
if lumped
% Gauss
[w, puntos, n] = gauss([3,3]);

Area = 0;
for i = 1:n
        Neta = [0, 1, 0 2*puntos(i,1), puntos(i,2), 0, 2*puntos(i,1)*puntos(i,2), puntos(i,2)^2]*A;%derivada de N en eta en los puntos de Gauss
        Nzeta = [0, 0, 1, 0, puntos(i,1), 2*puntos(i,2), puntos(i,1)^2, 2*puntos(i,1)*puntos(i,1)]*A;%derivada de N en zeta

        D = [Neta; Nzeta];

        J = D*nodos;
  
        mult = abs(det(J))*w(i);
        
        Area = Area + mult;
end% i
    m = Area*p/12*eye(cant_puntos);
    M(1:8,1:8) = m;%esquinas
    M(9:16,9:16) = 2*m;%lados
    return
end


%% Gauss
[w, puntos, n] = gauss([3,3]);

M=0;
dir1 = 1:2:8;
dir2 = 2:2:8;
for i = 1:n
        N = [1, puntos(i,1), puntos(i,2), puntos(i,1)^2, puntos(i,1)*puntos(i,2),  puntos(i,2)^2,  puntos(i,2)^2* puntos(i,2),  puntos(i,1)* puntos(i,2)^2]*A;

        Neta = [0, 1, 0 2*puntos(i,1), puntos(i,2), 0, 2*puntos(i,1)*puntos(i,2), puntos(i,2)^2]*A;%derivada de N en eta en los puntos de Gauss
        Nzeta = [0, 0, 1, 0, puntos(i,1), 2*puntos(i,2), puntos(i,1)^2, 2*puntos(i,1)*puntos(i,2)]*A;%derivada de N en zeta

        D = [Neta; Nzeta];

        J = D*nodos;
    
        Ns(dir1,1) = N;
        Ns(dir2,2) = N;
  
        mult = abs(det(J))*w(i);
        Mmin = Ns*p*Ns';
        
        M = M + Mmin*mult;
end% i
end