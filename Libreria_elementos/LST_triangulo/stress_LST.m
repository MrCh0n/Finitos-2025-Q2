function [Stress] = stress_LST(Coord, Uel, C, Czz)
%Stress = stress_LST(Coord, Uel, C, Czz)
%
%Devuelve las tension cte del elemento
%Devuelve un vector de 6 numeros
%
%Stress(1) es la tension xx
%Stress(2) es la tension yy
%Stress(3) es la tension xy
%Stress(4) una tension principal
%Stress(5) una tension principal
%Stress(6) una tension principal y zz
%
%las tensiones principales no estan ordenadas de mayor a menor
%
%Coord es una lista (x, y) de los 6 nodos del elemento
%
%Uel son los dezplamientos de los nodos, siendo un vector
%(Ux1,Uy1,Ux2,Uy2..) con el orden visto mas abajo
%
%C es el constitutivo de 3x3 y Czz la parte constitutiva de zz de 3x1
%
%poner los nodos como:
%    3
%   / \
%  6   5
% /     \
%1 - 4 - 2 
    Stress = zeros(1,6); %[xx yy xy principales] Stress(6) es zz
    
    % isoparametrico para tensiones LST
    cant_puntos = 6;
    
    x1 = [0; 1; 0; 0.5; 0; 0.5];
    y1 = [0; 0; 1; 0; 0.5; 0.5];
    A = [ones(cant_puntos,1) x1 y1 x1.*y1 x1.^2 y1.^2];
    
    %iso para extrapolar las tensiones en puntos de Gauss/superonvergentes
    x2 = [0; 1; 0];
    y2 = [0; 0; 1];
    A_rs = [ones(3,1) x2 y2];
    
    r = x1*2-1/3;%las esquinas de xi eta
    s = y1*2-1/3;

    %donde quiero las tensiones
    puntos =     [1/6 2/3
    1/6    1/6
    2/3    1/6];
    %toda la seccion superior se cambia junta

    dir1 = 1:2:2*cant_puntos;
    dir2 = 2:2:2*cant_puntos;
    
    orden = size(puntos,1);
    sigma_rs = zeros(orden,4);% (:,1) xx, 2 yy, 3 xy, 4 zz

    %tensiones en los puntos de Gauss
    for i = 1:orden
        Neta = [0, 1, 0, puntos(i,2), 2*puntos(i,1), 0]/A;%derivada de N en eta en los puntos de Gauss
        Nzeta = [0, 0, 1, puntos(i,1), 0, 2*puntos(i,2)]/A;%derivada de N en zeta
            
        D = [Neta; Nzeta];
    
        J = D*Coord;
    
        Bs = J\D;
            
        Bx = Bs(1,:);
        By = Bs(2,:);
    
        %crear la matrz B
        Bel(1,dir1) = Bx;
        Bel(2,dir2) = By;
        Bel(3,dir1) = By;
        Bel(3,dir2) = Bx;
       
        sigma_rs(i,1:3) = C*Bel*Uel; %sxx syy sxy
        sigma_rs(i,4) = Czz*Bel*Uel; %szz
    end

    sigma = zeros(1,4);
    for j = 1:cant_puntos
        N = [1 r(j) s(j)]/A_rs;
        sigma = sigma + N*sigma_rs;%sumo todos las tensiones en los nodos del elemento
    end

    sigma = sigma/cant_puntos;%se promedia

    Stress([1 2 3 6]) = sigma;

    sxx = Stress(1);
    syy = Stress(2);
    sxy = Stress(3);

    sigma_plano = [sxx sxy;
                   sxy syy];

    Stress(4:5) = eig(sigma_plano);
end