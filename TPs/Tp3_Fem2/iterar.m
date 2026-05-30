function [U, P, tiempo, presion] = iterar(U, P, K, Cg, F, Kp, M, M_monio, dt, freeU, freeP, R_est, nt)
    flag_relax = true;
    omega = 0.8;

    Kr = K(freeU,freeU);
    inv_Kr = inv(Kr);

    matriz_P = M/dt + Kp;
    
    if flag_relax
        matriz_P = matriz_P +M_monio/dt;
    end
    
    matriz_Pr = matriz_P(freeP,freeP);
    inv_Pr = inv(matriz_Pr);
    
    e = 1e-3;% 0.1%

    t = 0;
    presion = zeros(1,nt);
    tiempo = zeros(1,nt);

    for i = 1:nt
        presion(i) = P(1);
        tiempo(i) = t;
        fprintf("\n ---- Tiempo %.2f s ----\n",t);
        t = t+dt;

        difP = 1;
        U_old = U;
        P_old = P;
        n=0;
        while difP > e% Porciento de error
            %% Calculo U
            R = R_est + Cg*P;
    
            Rr = R(freeU);

            U_iter = U;
    
            U(freeU) = inv_Kr*Rr;%si haces un for loop ya esta hecho la inversion de matriz
            
            if flag_relax
                U = U*omega + U_iter*(1-omega);
            end


            %% Calculo P
            P_iter = P;
    
            Rp = M/dt*P_old + F*(U-U_old);%"fuerza en P"
            
            if flag_relax
                Rp = Rp + M_monio/dt*P_iter;
            end

            Rpr = Rp(freeP);%reducido
        
            P(freeP) = inv_Pr*Rpr;
            
            difP = abs(norm(P-P_iter)/norm(P));
            fprintf("%d %e\n",n, difP);
            n = n+1;
            if n > 300
                break
            end
        end

    end
end