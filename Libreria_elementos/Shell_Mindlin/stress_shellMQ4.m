function [Esfuerzos] = stress_shellMQ4(nodos, Uel, E,v,t)
%Stress = stress_shellMQ4(Coord, Uel, E,v,t)
%
%Devuelve las tensiones en el punto central
%Devuelve un vector de 6 numeros
%
%
%Momento(1) es Nx
%Momento(2) es Ny
%Momento(3) es Mx
%Momento(4) es My
%Momento(5) es Mxy
%Momento(6) es Qx
%Momento(7) es Qy
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
    L = zeros(5,6);
    L(1:3,1:3) = R;
    L(4:6,4:6) = [-v2;v1;v3];

    T = zeros(6*4);
    for i=1:4
        T((i-1)*6+1:i*6, (i-1)*6+1:i*6) = L;
    end
    Uelp = T*Uel;
    dir_placas = [3 4 5 9 10 11 15 16 17 21 22 23]; %dofs de placa
    dir_plano = [1 2 7 8 13 14 19 20]; %dofs de plano

    %Placa de Mindlin
    [Stress,Momentos] = stress_MQ4(nodosp(:,1:2), Uelp(dir_placas), E,v,t);

    C = E*t/(1-v^2)*[1 v 0;
                     v 1 0;
                     0 0 (1-v)/2];
    %Elemento membranal
    tension = stress_Q4(nodosp(:,1:2), Uelp(dir_plano), C, zeros(1,3));
    
    
    Nx = tension(1)*t; %sxx*t
    Ny = tension(2)*t; %syy*t

    Esfuerzos = [Nx Ny Momentos];
end