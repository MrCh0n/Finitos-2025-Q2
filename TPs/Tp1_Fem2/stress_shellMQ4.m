function [Stress,Momentos] = stress_shellMQ4(nodos, Uel, E,v,t)
%Stress = stress_shellMQ4(Coord, Uel, E,v,t)
%
%Devuelve las tensiones en el punto central
%Devuelve un vector de 6 numeros
%
%Stress(1) es la tension xx
%Stress(2) es la tension yy
%Stress(3) es la tension xy
%Stress(4) es la tension xz
%Stress(5) es la tensioon zy
%
%Momento(1) es Mx
%Momento(2) es My
%Momento(3) es Mxy
%Momento(4) es Qx
%Momento(5) es Qy
%
%
%Coord es una lista (x, y) de los 4 nodos del elemento
%
%Uel son los dezplamientos de los nodos, siendo un vector
%(w1,titax1,titay1,w2,titax2,...) con el orden visto mas abajo
%
%E, v y t son el modulo de elasticidad, coef de poisson, espesor
%
%poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 
    %% Matriz de rotacion
    v1 = (nodos(2,:)-nodos(1,:))/norm(nodos(2,:)-nodos(1,:));

    v3 = cross(v1,(nodos(4,:)-nodos(1,:)));

    v3 = v3 / norm(v3); % Normalize the cross product vector

    v2 = cross(v3,v1);

    R = [v1;v2;v3]; %v1 es vector fila

    %% Rotar nodos
    nodosp = (R*nodos')';

    %Placa de Mindlin
    Momento = stress_MQ4(nodosp(:,1:2), Uel, E,v,t);

    C = E*t/(1-v^2)*[1 v 0;
                     v 1 0;
                     0 0 (1-v)/2];
    %Elemento membranal
    tension = stress_Q4(nodosp(:,1:2), Uel, C, zeros(1,3));

    Nx = tension(1)*t; %sxx*t
    My = Momento(2);
    Qy = Momento(5);

end