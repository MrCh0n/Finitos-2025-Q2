function [K] = crearK_viga(nodos, E, I)
V = nodos(2,:) - nodos(1,:);
L = norm(V);

A = E*I/L;
Y1 = A*12/(L^2);
Y2 = A*6/L;
Y3 = A*4;
Y4 = A * 2; 

K =[Y1 Y2 -Y1 Y2;
    Y2 Y3 -Y2 Y4;
    -Y1 -Y2 Y1 -Y2;
    Y2 Y4 -Y2 Y3];
end