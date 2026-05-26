function [stress_el] = nodos_a_elem_Q4(stress)
% Promedia stress en los puntos de gauss para crear un stress promedio de
% elemento
    %% creo el isoparametrico
    cant_puntos = 4;
    
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(cant_puntos,1) x1 y1 x1.*y1];
    A=inv(A);
    %% Gauss
    [w, puntos, n] = gauss([2,2]);
    
    stress_el=0;
    for i = 1:n
            N = [1, puntos(i,1), puntos(i,2), puntos(i,1)*puntos(i,2)]*A;
      
            mult = w(i);
            Mmin = N*stress;
            
            stress_el = stress_el + Mmin*mult;
    end% i
    stress_el = stress_el/n;
end