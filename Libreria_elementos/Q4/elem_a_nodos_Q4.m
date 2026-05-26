function [stress_el] = elem_a_nodos_Q4(nodos, stress)
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
    
            % Neta = [0, 1, 0 puntos(i,2)]*A;
            % Nzeta = [0, 0, 1, puntos(i,1)]*A;
            % 
            % D = [Neta; Nzeta];
            % 
            % J = D*nodos;
      
            mult = w(i);
            Mmin = N'*stress;
            
            stress_el = stress_el + Mmin*mult;
    end% i
end