function [R] = carga_MQ4(nodos, carga, options)
%Crea el vector de cargas (Pz,Mx,My) de un elemento MQ4
%
%R = carga_MQ4(nodos, cargas_s, cargas_v)
%
%nodos:
%   es una lista (4,2) (x, y) de los 4 nodos del elemento
%
%cargas_s:
%   es la presion en z en cada uno de los nodos del MQ4
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
        carga (4,1) {mustBeNumeric}
    end



    %% creo el isoparametrico
    cant_puntos = 4;
    
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(cant_puntos,1) x1 y1 x1.*y1];
    A = inv(A);
    
    R = zeros(8,1);

    %para sacar dx/ds, dy/ds para cada linea
     tipo = [-1 0;
              0 -1;
              1 0;
              0 1;];%direccion de tau en la linea respectiva

     %giro la carga para que sea en x y
     if options.coordenada

        puntos_zeta = [y1(1:2)';
                        y1(2:3)';
                        y1([3 4])';
                        y1([4 1])'];
        
        puntos_eta = [x1(1:2)';
                        x1(2:3)';
                        x1(3:4)';
                        x1([4 1])'];

        posicion = [1 2;
                    3 4];

        for i = 1:4
            carga_s(i,:) = girar_carga(nodos, carga_s(i,:), A, puntos_eta(i,:), puntos_zeta(i,:), tipo(i,:), posicion);
        end
     end
    %% Gauss
    puntos = [-sqrt(3/5) 0 sqrt(3/5)];
    w = [5/9 8/9 5/9];
    
    orden = size(puntos,2);
    
    puntos_zeta = [-ones(orden,1)';
                    puntos;
                    ones(orden,1)';
                    puntos];
    
    puntos_eta = [puntos;
                   ones(orden,1)';
                   puntos
                   -ones(orden,1)'];
    
    nodos_linea = [1 2;
                   2 3;
                   3 4;
                   4 1];
    
    for i = 1:4
        tau = carga_s(i,1:2:end)';
        sigma = carga_s(i,2:2:end)';
        [Rx, Ry] = superficie(nodos, A, puntos_eta(i,:), puntos_zeta(i,:), nodos_linea(i,:), tau,sigma,w,tipo(i,:));
        R(2*nodos_linea(i,:) - 1) = R(2*nodos_linea(i,:) - 1) + Rx;
        R(2*nodos_linea(i,:)) = R(2*nodos_linea(i,:)) + Ry;
    end
    
    volumen =   area(nodos,puntos,w,A); %es de 4x4
    %volumen*[1;1;1;1;1;1;1;1] la fuerza/m^2 en cada nodo, normalmente rho*g*t
    
    R(1:2:end) = R(1:2:end) + volumen*carga_v(:,1);
    R(2:2:end) = R(2:2:end) + volumen*carga_v(:,2);
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


