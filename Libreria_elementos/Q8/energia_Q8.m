function [U] = energia_Q8(nodos,epsilon, C)
%Devuelve la energia de deformacion U de un Q8
%
%U = energia_Q8(nodos, epsilon, C)
%
%nodos es una lista (x, y) de los 8 nodos del elemento
%
%epsilon es un vector [exx1, eyy1, exy1, exx2..., exy8] como el de
%deformaciones_Q8
%
%C es el constitutivo de 3x3
%
%poner los nodos como:
% 4 - 7 - 3
% |       |
% 8       6
% |       |
% 1 - 5 - 2 
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

%% Gauss e*C*e
exx = epsilon(1:3:end)';
eyy = epsilon(2:3:end)';
exy = epsilon(3:3:end)';

U = 0;

[w, puntos, n] = gauss([3,3]);

for i = 1:n
    xi = puntos(i,1);
    eta = puntos(i,2);
    
    N = [1, xi, eta, xi^2, xi*eta,  eta^2,  xi^2*eta,  xi*eta^2]*A;

    Neta = [0, 1, 0, 2*xi, eta, 0, 2*xi*eta, eta^2]*A;%derivada de N en eta en los puntos de Gauss
    Nzeta = [0, 0, 1, 0, xi, 2*eta, xi^2, 2*xi*eta]*A;%derivada de N en zeta

    D = [Neta; Nzeta];

    J = D*nodos;

    e_xx = N*exx;
    e_yy = N*eyy;
    e_xy = N*exy;

    e = [e_xx; e_yy; e_xy];

    mult = abs(det(J))*w(i);
    integrando = e'*C*e;
    
    U = U + integrando*mult;
end% i

end