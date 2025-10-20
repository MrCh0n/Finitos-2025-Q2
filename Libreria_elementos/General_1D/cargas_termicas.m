function [sigma_t, Fuerzas] = cargas_termicas(A,E,dt,alpha)
    sigma_t = -alpha*E*dt;

    Fuerzas = sigma_t*[1; -1]*A;
end