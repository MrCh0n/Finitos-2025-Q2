function [Uel, Vel] = flujo_Q4(coord, Tel, k)
    % isoparametrico
    cant_puntos = 4;
    
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(cant_puntos,1) x1 y1 x1.*y1];
    
    %iso para extrapolar las tensiones en puntos de Gauss/superonvergentes
    x2 = [-1; 1; 1; -1];
    y2 = [-1; -1; 1; 1];
    A_rs = [ones(4,1) x2 y2 x2.*y2];
    
    r = x1*sqrt(3);%las esquinas de xi eta
    s = y1*sqrt(3);

    %donde quiero la tensiones
    puntos = [-sqrt(1/3) sqrt(1/3)];
     %toda la seccion superior se cambia junta
    
    cant = size(puntos,2);
    flujo_rs = zeros(cant^2,2);% (:,1) xx, 2 yy, 3 xy, 4 zz
    
    %tension en los puntos de Gauss
    for i = 1:cant
        offset = (i-1)*cant;
        for j = 1:cant
            Neta = [0, 1, 0 puntos(j)]/A;
            Nzeta = [0, 0, 1, puntos(i)]/A;

            D = [Neta; Nzeta];

            J = D*coord;

            Bs = J\D;
        
            flujo_rs(offset+j,:) = -k*Bs*Tel; %(qx,qy)
        end
    end
    flujo_rs = flujo_rs([1 2 4 3],:);%queda mal ubicado del for loop puntos(i/j)
   

    Uel = 0;
    Vel = 0;
    
    for j = 1:cant_puntos
        N = [1 r(j) s(j) r(j)*s(j)]/A_rs;
        Uel = Uel + N*flujo_rs(:,1);%sumo todos los flujos en x en los nodos del elemento
        Vel = Vel + N*flujo_rs(:,2);%sumo todos los flujos en y los nodos del elemento
    end

    Uel = Uel/cant_puntos;
    Vel = Vel/cant_puntos;
end