function [K] = crearK_viga_barra(nodos,E,A,I)
%crea la K local a los nodos dados
%
%K = crearK_viga_barra(nodos, E, A, I)
%
%nodos es una lista (x, y) de los 2 nodos del elemento
%E es el modulo de Young del materia
%A es el area de la seccion
%I es el momento de Inercia de la seccion

    V = nodos(2,:) - nodos(1,:);
    L = norm(V);

    Kel_bar = A*E/L *[1 -1;-1 1];
    
    Y1 = 12*E*I/L^3;
    Y2 = 6*E*I/L^2;
    Y3 = 4*E*I/L;
    Y4 = 2*E*I/L;

    Kel_vig = [Y1 Y2 -Y1 Y2;
           Y2 Y3 -Y2 Y4;
           -Y1 -Y2 Y1 -Y2;
           Y2 Y4 -Y2 Y3];

    dir_bar = [1 4];
    dir_vig = [2 3 5 6];

    K(dir_bar, dir_bar) = Kel_bar;
    K(dir_vig, dir_vig) = Kel_vig;

    cs = V/L;
    c = cs(1);
    s = cs(2);

    Q = [c s 0;
         -s c 0;
         0 0 1];
    T(1:3, 1:3) = Q;
    T(4:6, 4:6) = Q;
    
    K = T'*K*T;

end