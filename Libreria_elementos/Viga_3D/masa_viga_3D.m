function [M] = masa_viga_3D(nodos, rho, A, Ip, auxiliar)
    V = nodos(2,:) - nodos(1,:);
    L = norm(V);

    Mel([1 7],[1 7]) = rho*A*L/2*[1 0; 0 1];
    Mel([2 6 8 12],[2 6 8 12]) = rho*A*L/24 * diag([12 L^2 12 L^2]);
    Mel([3 5 9 11],[3 5 9 11]) = rho*A*L/24 * diag([12 L^2 12 L^2]);
    Mel([4 10],[4 10]) = rho*Ip*L/2*[1 0; 0 1];

    dir1 = V/L;
    Q = crearQ(dir1, auxiliar);

    M = Q'*Mel*Q;
end

function [Q] = crearQ(dir1, auxiliar)
    dir3 = cross(auxiliar, dir1)/norm(cross(auxiliar, dir1));
    dir2 = cross(dir3, dir1);
    lambda = [dir1; dir2; dir3];
    assert(abs(det(lambda))-1 < 1e-6);
    
    Q = blkdiag(lambda, lambda, lambda, lambda);
end