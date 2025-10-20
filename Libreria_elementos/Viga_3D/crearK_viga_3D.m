function [K] = crearK_viga_3D(nodos,E,G,A,Iy,Iz,Ip,auxiliar)
    V = nodos(2,:) - nodos(1,:);
    L = norm(V);
    
    Kel = zeros(12);
    Kel([1 7],[1 7]) = A*E/L*[1 -1;-1 1];
    Kel([2 6 8 12],[2 6 8 12]) = rigidez_viga(L,E,Iz);
    Kel([3 5 9 11],[3 5 9 11]) = rigidez_viga(L,E,Iy);
    Kel([4 10],[4 10]) = G*Ip/L*[1 -1;-1 1];

    dir1 = V/L;
    Q = crearQ(dir1, auxiliar);
    K = Q'*Kel*Q;
end

function [Q] = crearQ(dir1, auxiliar)
    dir3 = cross(auxiliar, dir1)/norm(cross(auxiliar, dir1));
    dir2 = cross(dir3, dir1);
    lambda = [dir1; dir2; dir3];
    assert(abs(det(lambda))-1 < 1e-6);
    
    Q = blkdiag(lambda, lambda, lambda, lambda);
end

function [K] = rigidez_viga(L, E, I)
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