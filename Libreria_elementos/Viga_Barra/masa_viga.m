function [M] = masa_viga(rho, Area, L, lumped)
    syms x
    A_barra = [1 0;
               1 L];

    A_viga = [1 0 0 0;
              0 1 0 0;
              1 L L^2 L^3;
              0 1 2*L 3*L^2];

    M = zeros(6);
    dir_barra = [1 4];
    dir_viga = [2 3 5 6];

    %Gauss
    cant_puntos = 4;

    a = 0.339981;%ver cuadratura de 4
    b = 0.861136;
    puntos = [-a a -b b]*L/2 + L/2;

    wa = 0.652145;
    wb = 1-wa;
    w = [wa wa wb wb]*L/2;
    area = 0;
    for i = 1:cant_puntos
        N_barra = [1 puntos(i)]/A_barra;
        N_viga = [1 puntos(i) puntos(i)^2 puntos(i)^3]/A_viga;

        M(dir_viga,dir_viga) = M(dir_viga,dir_viga) + N_viga'*N_viga*w(i);
        M(dir_barra,dir_barra) = M(dir_barra,dir_barra) + N_barra'*N_barra*w(i);

    end
    M = M*rho;

    if lumped%no funciona sumar las filas para la viga, para la barra si
        aux = zeros(6);
        for i = 1:6
            aux(i,i) = sum(M(i,:));
        end
        M = aux;
    end
end