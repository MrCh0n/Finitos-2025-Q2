function [Ks] = pandeo_Ks_viga(L)
    A = [1 0 0 0;
        0 1 0 0;
        1 L L^2 L^3;
        0 1 2*L 3*L^2];

    Ks = zeros(6);
    dir_viga = [2 3 5 6];
    %Gauss
    cant_puntos = 3;
  
    puntos = [-sqrt(3/5) 0 sqrt(3/5)]*L/2 + L/2;
    w = [5/9 8/9 5/9]*L/2;
    for i = 1:cant_puntos
        G = [0 1 2*puntos(i) 3*puntos(i)^2]/A;
        Ks(dir_viga,dir_viga) = Ks(dir_viga,dir_viga) + G'*G*w(i);
    end
end