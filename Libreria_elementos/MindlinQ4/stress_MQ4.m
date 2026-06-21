function [Stress,Momentos] = stress_MQ4(Coord, Uel, E,v,t)
%Stress = stress_Q4(Coord, Uel, C, Czz)
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

    %% matriz constitutiva
    G = E/(2*(1+v));
    D = E*t^3/(12*(1-v^2));
    k = 5/6;
    
    D_b = [D v*D 0;
          v*D D 0;
          0 0 (1-v)/2*D];
    D_s = [k*G*t 0;
          0 k*G*t];

    C_b = E/(1-v^2)*[1 v 0;
                     v 1 0;
                     0 0 (1-v)/2];
    C_s = G*[1 0; 0 1];

    %% Isoparametrico
    cant_puntos = 4;
    
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(cant_puntos,1) x1 y1 x1.*y1];

    B_b = zeros(3,3*cant_puntos);
    B_s = zeros(2,3*cant_puntos);
    
    dir1 = 1:3:3*cant_puntos;
    dir2 = 2:3:3*cant_puntos;
    dir3 = 3:3:3*cant_puntos;
    
    % uso solo un punto, el (0,0) para evaluar B y tensiones
    N = [1,0,0,0]/A;
    Neta = [0 1 0 0]/A;
    Nzeta = [0 0 1 0]/A;
    
    D = [Neta; Nzeta];

    J = D*Coord;

    Bs = J\D;
        
    Bx = Bs(1,:);
    By = Bs(2,:);

    %crear la matrz B
    B_b(1,dir2) = Bx;
    B_b(2,dir3) = By;
    B_b(3,dir2) = By;
    B_b(3,dir3) = Bx;

    B_s(1,dir1) = -Bx;
    B_s(2,dir1) = -By;
    B_s(1,dir2) = N;
    B_s(2,dir3) = N;

    %% Tension
    Stress = zeros(1,5); 

    z = t/2; % calculamos a altura t/2 (deberiamos tmb calcularlo a -t/2)
    % OOJJJOOOOO

    eps_b = B_b*Uel;
    eps_s = B_s*Uel;
    eps = [eps_b;eps_s];

    S = diag([-z -z -z 1 1]);
    C = zeros(5);
    C(1:3,1:3) = C_b;
    C(4:5,4:5)= C_s;

    Stress = t*C*S*eps; %Nx Ny xy

    Momentos(1:3) = D_b*eps(1:3);

    Momentos(4:5) = D_s*eps(4:5);
end