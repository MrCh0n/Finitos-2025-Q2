function [M] = masa_viga(rho, Area, L, lumped)
    M = zeros(6);
    dir_barra = [1 4];
    dir_viga = [2 3 5 6];

    masa = rho*Area*L;
    
    if lumped
        M_barra = [1 0;
                   0 1]/2;

        M_viga = [12 0 0 0;
                  0 L^2 0 0;
                  0 0 12 0;
                  0 0 0 L^2]/24;
    else
        M_barra = [2 1;
                   1 2]/6;
    
        M_viga = [156 22*L 54 -13*L;
                  22*L 4*L^2 13*L -3*L^2;
                  54 13*L 156 -22*L;
                  -13*L -3*L^2 -22*L 4*L^2];
    end

    M(dir_barra,dir_barra) = M_barra;
    M(dir_viga, dir_viga) = M_viga;

    M = M*masa;
end