function [stress] = stress_viga_3d(nodos, U ,E,G,y_max,z_max,r,auxiliar, div)
    V = nodos(2,:) - nodos(1,:);
    L = norm(V);

    dir1 = V/L;
    Q = crearQ(dir1, auxiliar);

    U_local = Q*U;

    %flexion
    x = linspace(0,L,div);
    Bflexion = [-6/L^2+12*x/L^3
                -4/L+6*x/L^2
                 6/L^2-12*x/L^3
                -2/L+6*x/L^2];

    stressAx = E*(U_local(7)-U_local(1))/L;

    strain1 = Bflexion'*U_local([2 6 8 12]);   %1 es con v1, titaz1, v2, titaz2
    stress1 = E*y_max*strain1;

    strain2 = Bflexion'*U_local([3 5 9 11]);   %2 es con w1, titay1, w2, titay2
    stress2 = E*z_max*strain2;

    stressFlex = max(sqrt(stress1.^2 + stress2.^2));

    stressxx = max(abs([stressFlex+stressAx stressFlex-stressAx]));

    stressxy = G*r*(U_local(10)-U_local(4))/L;
    
    vonMinses = sqrt(stressxx^2 + 3*stressxy^2);

    stress(1:5) = [vonMinses, stressAx, stressFlex,stressxx, stressxy]/1e6;%los paso a Mpa
end

function [Q] = crearQ(dir1, auxiliar)
    dir3 = cross(auxiliar, dir1)/norm(cross(auxiliar, dir1));
    dir2 = cross(dir3, dir1);
    lambda = [dir1; dir2; dir3];
    assert(abs(det(lambda))-1 < 1e-6);
    
    Q = blkdiag(lambda, lambda, lambda, lambda);
end