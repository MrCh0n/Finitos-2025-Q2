function [M] = masa_viga(nodos, rho, A)
% lumpeada
    V = nodos(2,:) - nodos(1,:);
    L = norm(V);

    M = rho*A*L/24 * diag([12 L^2 12 L^2]);
end