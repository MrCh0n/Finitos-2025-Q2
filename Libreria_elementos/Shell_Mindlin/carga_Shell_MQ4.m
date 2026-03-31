function [R] = carga_Shell_MQ4(nodos, carga_v) %TODO
%Crea el vector de cargas (Pz,Mx,My) de un elemento Shell MQ4
%
%R = carga_MQ4(nodos, cargas_s, cargas_v)
%
%nodos:
%   es una lista (4,2) (x, y) de los 4 nodos del elemento
%
%cargas_v:
%   es la presion/esfuerzo en z en cada uno de los nodos del shell_MQ4, tmb puede
%   ser volumetrica (premultiplicar por t)
%
%
%
% poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 

arguments
        nodos (4,2) {mustBeNumeric}
        carga_v (4,1) {mustBeNumeric}
    end
    
    % isoparametrico
    cant_puntos = 4;
    
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(cant_puntos,1) x1 y1 x1.*y1];
    A = inv(A);

    R = zeros(12,1);

    % Gauss
    puntos = [-sqrt(3/5) 0 sqrt(3/5)];
    w = [5/9 8/9 5/9];
    
    volumen =   area(nodos,puntos,w,A); %es de 4x4
    %volumen*[1;1;1;1;1;1;1;1] la fuerza/m^2 en cada nodo, normalmente rho*g*t
    
    R(3:6:end) = volumen*carga_v(:,1);
end



function [F_volumen] = area(nodos,puntos,w,A)
% Gauss
orden = size(puntos,2);

F_volumen = zeros(4);

for i = 1:orden
    for j = 1:orden
        N = [1 puntos(i), puntos(j), puntos(i)*puntos(j)]*A;%N
        Neta = [0, 1, 0, puntos(j)]*A;%derivada de N en eta en los puntos de Gauss
        Nzeta = [0, 0, 1, puntos(i)]*A;%derivada de N en zeta
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
        
        integrando = N'*N;
    
        F_volumen = F_volumen + integrando*det(J)*w(i)*w(j);
    end% j
end% i

end
