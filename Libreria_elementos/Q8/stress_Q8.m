function [Stress] = stress_Q8(Coord, Uel, C, Czz)
%Stress = stress_Q8(Coord, Uel, C, Czz)
%
%Devuelve las tensiones promediadas de los nodos
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
%Coord es una lista (x, y) de los 8 nodos del elemento
%
%Uel son los dezplamientos de los nodos, siendo un vector
%(Ux1,Uy1,Ux2,Uy2,Ux3...) con el orden visto mas abajo
%
%C es el constitutivo de 3x3 y Czz la parte constitutiva de zz de 3x1
%
%poner los nodos como:
% 4 - 7 - 3
% |       |
% 8       6
% |       |
% 1 - 5 - 2 

    Stress = zeros(1,6); %[xx yy xy principales] Stress(6) es zz
    
    % isoparametrico
    cant_puntos = 8;
    
    x1 = [-1; 1; 1; -1; 0; 1; 0; -1];
    y1 = [-1; -1; 1; 1; -1; 0; 1; 0];
    A = [ones(cant_puntos,1) x1 y1 x1.^2 x1.*y1 y1.^2 x1.^2.*y1 y1.^2.*x1];
    
    %iso para extrapolar las tensiones en puntos de Gauss/superonvergentes
    x2 = [-1; 1; 1; -1];
    y2 = [-1; -1; 1; 1];
    A_rs = [ones(4,1) x2 y2 x2.*y2];
    
    r = x1*sqrt(3);%las esquinas de xi eta
    s = y1*sqrt(3);

    %donde quiero la tensiones
    puntos = [-sqrt(1/3) sqrt(1/3)];
     %toda la seccion superior se cambia junta

    dir1 = 1:2:2*cant_puntos;
    dir2 = 2:2:2*cant_puntos;
    
    cant = size(puntos,2);
    sigma_rs = zeros(cant^2,4);% (:,1) xx, 2 yy, 3 xy, 4 zz
    
    %tension en los puntos de Gauss
    for i = 1:cant
        offset = (i-1)*cant;
        for j = 1:cant
            Neta = [0, 1, 0 2*puntos(i), puntos(j), 0, 2*puntos(i)*puntos(j), puntos(j)^2]/A;%derivada de N en eta en los puntos de Gauss
            Nzeta = [0, 0, 1, 0, puntos(i), 2*puntos(j), puntos(i)^2, 2*puntos(i)*puntos(j)]/A;%derivada de N en zeta
                
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
        
            sigma_rs(offset+j,1:3) = C*Bel*Uel; %sxx syy sxy
            sigma_rs(offset+j,4) = Czz*Bel*Uel; %szz
        end
    end
    
    sigma = zeros(1,4);
    for j = 1:cant_puntos
        N = [1 r(j) s(j) r(j)*s(j)]/A_rs;
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