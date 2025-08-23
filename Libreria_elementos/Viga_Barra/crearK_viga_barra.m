function [K] = crearK_viga_barra(nodos,E,A,I)
%crea la K local a los nodos dados
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

end