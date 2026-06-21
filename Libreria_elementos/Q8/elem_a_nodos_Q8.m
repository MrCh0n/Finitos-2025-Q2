function [stress_el] = elem_a_nodos_Q8(nodos, stress)
    %% creo el isoparametrico
    cant_puntos = 8;
    
    puntos_Q8 = [-1    -1;
                  1    -1;
                  1     1;
                 -1     1;
                  0    -1;
                  1     0;
                  0     1;
                 -1     0];
    
    x1 = puntos_Q8(:,1);%puntos del cuadrado
    y1 = puntos_Q8(:,2);
    A = [ones(cant_puntos,1) x1 y1 x1.^2 x1.*y1 y1.^2 x1.^2.*y1 y1.^2.*x1];
    A = inv(A);
    %% Gauss
    [w, puntos, n] = gauss([2,2]);
    
    stress_el=0;
    for i = 1:n
        xi = puntos(i,1);
        eta = puntos(i,2);
    
        N = [1, xi, eta, xi^2, xi*eta,  eta^2,  xi^2*eta,  xi*eta^2]*A;
    
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