function P_ana = presion_terzaghi(y, t, P0, H, cv)
    % Evalúa la serie analítica de Terzaghi para una altura y tiempo dados
    tmp = 0;
    if t == 0
        P_ana = P0;
        return;
    end
    
    cte = pi^2 * t * cv / (H^2 * 4);
    for k = 1:1000 % Truncamiento de la sumatoria
        exponencial = exp(-(2*k-1)^2 * cte);
        tmp = tmp + (-1)^(k-1)/(2*k-1) * cos(pi/2 * y/H * (2*k-1)) * exponencial;
    end
    P_ana = 4/pi * P0 * tmp;
end