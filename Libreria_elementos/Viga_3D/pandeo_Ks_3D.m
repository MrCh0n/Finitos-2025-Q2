function [Ks] = pandeo_Ks_3D(nodos, auxiliar)
    V = nodos(2,:) - nodos(1,:);
    L = norm(V);

    Ksel = zeros(12);
    Ksel([2 6 8 12],[2 6 8 12]) = 1/(30*L)*[36 3*L -36 3*L;
                                          3*L 4*L^2 -3*L -L^2;
                                          -36 -3*L 36 -3*L;
                                          3*L -L^2 -3*L 4*L^2];
    Ksel([3 5 9 11],[3 5 9 11]) = 1/(30*L)*[36 3*L -36 3*L;
                                          3*L 4*L^2 -3*L -L^2;
                                          -36 -3*L 36 -3*L;
                                          3*L -L^2 -3*L 4*L^2];

    dir1 = V/L;
    Q = crearQ(dir1, auxiliar);

    Ks = Q'*Ksel*Q;
end

function [Q] = crearQ(dir1, auxiliar)
    dir3 = cross(auxiliar, dir1)/norm(cross(auxiliar, dir1));
    dir2 = cross(dir3, dir1);
    lambda = [dir1; dir2; dir3];
    assert(abs(det(lambda))-1 < 1e-6);
    
    Q = blkdiag(lambda, lambda, lambda, lambda);
end