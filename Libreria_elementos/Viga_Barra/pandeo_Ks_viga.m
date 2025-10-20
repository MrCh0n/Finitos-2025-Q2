function [Ks] = pandeo_Ks_viga(nodos)
    V = nodos(2,:) - nodos(1,:);
    L = norm(V);

    Ks = zeros(6);
    dir_viga = [2 3 5 6];
    %Gauss
   
    Ks(dir_viga,dir_viga) = 1/(30*L)*[36 3*L -36 3*L;
                              3*L 4*L^2 -3*L -L^2;
                              -36 -3*L 36 -3*L;
                              3*L -L^2 -3*L 4*L^2];
    cs = V/L;
    c = cs(1);
    s = cs(2);

    Q = [c s 0;
         -s c 0;
         0 0 1];
    T(1:3, 1:3) = Q;
    T(4:6, 4:6) = Q;

    Ks = T'*Ks*T;
end