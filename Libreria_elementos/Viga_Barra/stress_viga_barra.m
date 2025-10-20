function [stress, von_Minses] = stress_viga_barra(nodos, U ,E ,y_max ,div, stress_termico)
    V = nodos(2,:) - nodos(1,:);
    L = norm(V);

    cs = V/L;
    c = cs(1);
    s = cs(2);

    Q = [c s 0;
         -s c 0;
         0 0 1];

    T(1:3, 1:3) = Q;
    T(4:6, 4:6) = Q;

    U_local = T*U;

    %flexion
    x = linspace(0,L,div);
    Bflexion = [-6/L^2+12*x/L^3
                -4/L+6*x/L^2
                 6/L^2-12*x/L^3
                -2/L+6*x/L^2];

    stressAx = E*(U_local(4)-U_local(1))/L + stress_termico;

    strain1 = Bflexion'*U_local([2 3 5 6]);   %1 es con v1, titaz1, v2, titaz2
    stress1 = E*y_max*strain1;

    stressFlex = max(sqrt(stress1.^2));

    stressxx = abs([stressFlex+stressAx stressFlex-stressAx]);
    
    stressxy = 0;
    
    von_Minses = abs(stress1)+abs(stressAx);

    stressxx = max(stressxx);

    vonMinses = sqrt(stressxx^2 + 3*stressxy^2);

    stress(1:5) = [vonMinses, stressAx, stressFlex,stressxx, stressxy]/1e6;%los paso a Mpa
end