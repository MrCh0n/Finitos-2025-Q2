function [Stress] = stress_Q4(Coord, Uel, C, Czz)
%Stress = stress_Q4(Coord, Uel, C, Czz)
%
%Devuelve las tensiones en el punto central
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
%Coord es una lista (x, y) de los 4 nodos del elemento
%
%Uel son los dezplamientos de los nodos, siendo un vector
%(Ux1,Uy1,Ux2,Uy2,Ux3...) con el orden visto mas abajo
%
%C es el constitutivo de 3x3 y Czz la parte constitutiva de zz de 3x1
%
%poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 

    Stress = zeros(1,6); %[xx yy xy principales] Stress(6) es zz
    
    % isoparametrico
    cant_puntos = 4;
    
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(cant_puntos,1) x1 y1 x1.*y1];
    Bel = zeros(3,2*cant_puntos);
        
    dir1 = 1:2:2*cant_puntos;
    dir2 = 2:2:2*cant_puntos;
    
    % uso solo un punto, el (0,0) para evaluar B y tensiones
    Neta = [0 1 0 0]/A;
    Nzeta = [0 0 1 0]/A;
    
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

    Stress(1:3) = C*Bel*Uel; %sxx syy sxy
    Stress(6) = Czz*Bel*Uel; %szz

    sxx = Stress(1);
    syy = Stress(2);
    sxy = Stress(3);

    sigma_plano = [sxx sxy;
                   sxy syy];

    Stress(4:5) = eig(sigma_plano);
end