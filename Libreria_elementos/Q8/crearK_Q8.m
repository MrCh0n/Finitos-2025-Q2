function [K] = crearK_Q8(nodos, C)
%Crea la matriz de resistencia de un Q8 sin espesor
%
%K = crearK_Q8(nodos, C)
%
%nodos es una lista (x, y) de los 8 nodos del elemento
%C es el constitutivo de 3x3
%
%poner los nodos como:
% 4 - 7 - 3
% |       |
% 8       6
% |       |
% 1 - 5 - 2 
arguments (Input)
   nodos(8,2) {mustBeNumeric}
   C(3,3) {mustBeNumeric}
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


%% Gauss
puntos = [-sqrt(3/5) 0 sqrt(3/5)];
w = [5/9 8/9 5/9];

orden = size(puntos,2);

B = zeros(3,2*cant_puntos);
K = 0;

for i = 1:orden
    for j = 1:orden
        Neta = [0, 1, 0 2*puntos(i), puntos(j), 0, 2*puntos(i)*puntos(j), puntos(j)^2]/A;%derivada de N en eta en los puntos de Gauss
        Nzeta = [0, 0, 1, 0, puntos(i), 2*puntos(j), puntos(i)^2, 2*puntos(i)*puntos(j)]/A;%derivada de N en zeta
        
        J = [Neta; Nzeta]*nodos;%el jacobiaos para el punto de Gauss

        J1 = inv(J);%el inverso para hacer B
    
        Bx = Neta*J1(1,1) + Nzeta*J1(1,2);%transforma las derivadas a x
        By = Neta*J1(2,1) + Nzeta*J1(2,2);%transformo las derivadas a y
        
        dir1 = 1:2:2*cant_puntos;
        dir2 = 2:2:2*cant_puntos;
        %crear la matrz B
        B(1,dir1) = Bx;
        B(2,dir2) = By;
        B(3,dir1) = By;
        B(3,dir2) = Bx;
    
        mult = abs(det(J))*w(i)*w(j);%pesos y cambio de area abs por si es negativo el Jacobiano
        Kmin = B'*C*B;%matriz a integrar
        
        K = K + Kmin*mult;

    end% j
end% i

end