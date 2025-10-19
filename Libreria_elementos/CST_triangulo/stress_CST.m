function [Stress] = stress_CST(Coord, Uel, C, Czz)
%Stress = stress_CST(Coord, Uel, C, Czz)
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
%Coord es una lista (x, y) de los 3 nodos del elemento
%
%Uel son los dezplamientos de los nodos, siendo un vector
%(Ux1,Uy1,Ux2,Uy2...) con el orden visto mas abajo
%
%C es el constitutivo de 3x3 y Czz la parte constitutiva de zz de 3x1
%
%poner los nodos como:
%    3
%   / \
%  /   \
% /     \
%1 - - - 2 
    Stress = zeros(1,6); %[xx yy xy principales] Stress(6) es zz
    
    x1 = Coord(1,1);
    y1 = Coord(1,2);
    x2 = Coord(2,1);
    y2 = Coord(2,2);
    x3 = Coord(3,1);
    y3 = Coord(3,2);
    
    A = det([1 x1 y1; 1 x2 y2;1 x3 y3]);

    Bel = 1/A*[y2-y3 0    y3-y1 0     y1-y2 0;
               0   x3-x2   0     x1-x3 0   x2-x1;
               x3-x2 y2-y3 x1-x3 y3-y1 x2-x1 y1-y2];

    Stress(1:3) = C*Bel*Uel; %sxx syy sxy
    Stress(6) = Czz*Bel*Uel; %szz

    sxx = Stress(1);
    syy = Stress(2);
    sxy = Stress(3);

    sigma_plano = [sxx sxy;
                   sxy syy];

    Stress(4:5) = eig(sigma_plano);
end