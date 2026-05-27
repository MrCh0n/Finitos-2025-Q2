function [U, P, tiempo, presion] = iterar(U, P, K, Cg, F, Kp, M, dt, freeU, freeP, R_est, iter)
    Kr = K(freeU,freeU);
    inv_Kr = inv(Kr);

    matriz_P = M/dt + Kp;
    matriz_Pr = matriz_P(freeP,freeP);
    inv_Pr = inv(matriz_Pr);
    
    e = 10e-2;% 5%

    t = 0;
    presion = zeros(1,iter);
    tiempo = zeros(1,iter);

    for i = 1:iter
        presion(i) = norm(P);
        tiempo(i) = t;
        t = t+dt;

        difP = 1;
        U_iter = U;
        P_iter = P;
        n=0;
        while difP > e% Porciento de error
            %% Calculo U
            R = R_est - Cg*P;
    
            Rr = R(freeU);

            %U_prev = U;
    
            U(freeU) = inv_Kr*Rr;%si haces un for loop ya esta hecho la inversion de matriz          
            %% Calculo P
            P_prev = P;
    
            Rp = M*P_iter + F*(U-U_iter);%"fuerza en P"
            Rpr = Rp(freeP);%reducido
    
            P(freeP) = inv_Pr*Rpr;
            
            difP = abs(norm(P-P_prev)/norm(P));
            fprintf("%d %f\n",n, difP);
            n = n+1;
        end

    end
end