function [U] = energia_Q4(nodos,epsilon, C)
%Devuelve la energia de deformacion U de un Q4
%
%U = energia_Q4(nodos, epsilon, C)
%
%nodos es una lista (x, y) de los 4 nodos del elemento
%
%epsilon es un vector [exx1, eyy1, exy1, exx2..., exy4] como el de
%deformaciones_Q4
%
%C es el constitutivo de 3x3
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
%% Gauss e*C*e
exx = epsilon(1:3:end)';
eyy = epsilon(2:3:end)';
exy = epsilon(3:3:end)';

U = 0;

[w, puntos, n] = gauss([2,2]);

for i = 1:n
    xi = puntos(i,1);
    eta = puntos(i,2);

    N = [1, xi, eta, xi*eta]*A;

    Neta = [0, 1, 0 eta]*A;
    Nzeta = [0, 0, 1, xi]*A;

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