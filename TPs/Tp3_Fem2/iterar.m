function [U, P] = iterar(U, P, K, Cg, F, Kp, M, freeU, freeP, R_est, iter)
    Kr = K(freeU,freeU);
    inv_Kr = inv(Kr);

    matriz_P = M/dt + Kp;
    matriz_Pr = matriz_P(freeP,freeP);
    inv_Pr = inv(matriz_Pr);
    
    for i = 1:iter
        %% Calculo U
        R = R_est - Cg*P;

        Rr = R(free);
        
        U_prev = U;

        U(free) = inv_Kr*Rr;%si haces un for loop ya esta hecho la inversion de matriz
        
        difU = norm(U-U_prev);
        
        %% Calculo P
        P_prev = P;

        Rp = M*P_prev + F*(U-U_prev);%"fuerza en P"
        Rpr = Rp(freeP);%reducido

        P(freeP) = inv_Pr*Rpr;
        
        difP = norm(P-P_prev);
    end
end