function [K] = crearK_shell_degeneradoQ4(nodos,E,v,t,options)
%Crea la matriz de resistencia de una cascara degenerada Q4
%
%K = crearK_shell_degenerado(nodos, E, v, t)
%
%nodos es una lista (x, y, z) de los 4 nodos del elemento
%E es el modulo de elasticidad
%v el coef de poisson
%t una lista (1,4) del espesor del nodo correspondiente
%
%poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 

arguments (Input)
   nodos(4,3) {mustBeNumeric}
   E {mustBeNumeric}
   v {mustBeNumeric}
   t {mustBeNumeric}
   options.V=[];
end
%cantidad de dofs por nodo
dofsxnod = 5;
%% Constitutivo de degenerado
D1 = E/(1 - v^2);
G = 0.5*E/(1 + v);
c = 5/6;
D = [   D1 v*D1  0   0   0
      v*D1   D1  0   0   0
          0     0  G   0   0
          0     0  0 c*G   0
          0     0  0   0 c*G ];

Ds = zeros(size(D));
Ds(1:3,1:3)= D(1:3,1:3);
Dc = zeros(size(D));
Dc(4:end,4:end) = D(4:end,4:end);
%% creo el isoparametrico
cant_puntos = 4;

x1 = [-1; 1; 1; -1];
y1 = [-1; -1; 1; 1];
A = [ones(cant_puntos,1) x1 y1 x1.*y1];
A = inv(A);

%% Direcciones de los nodos
v = zeros(3,3,cant_puntos);
for i = 1:cant_puntos
    switch i
        case 1
            x1=1; x2=2; y1=1; y2=4;
        case 2
            x1=1; x2=2; y1=2; y2=3;
        case 3
            x1=4; x2=3; y1=2; y2=3;
        case 4
            x1=4; x2=3; y1=1; y2=4;
    end
    v1 = (nodos(x2,:)-nodos(x1,:))/norm(nodos(x2,:)-nodos(x1,:));

    v3 = cross(v1,(nodos(y2,:)-nodos(y1,:)));

    v3 = v3 / norm(v3); % Normalize the cross product vector

    v2 = cross(v3,v1);

    v(:,:,i) = [v1;v2;v3]'; %v1 es vector fila
end%i
%% Gauss Flexion
[w, puntos, n] = gauss([2,2,2]);

%para hacer el jacobiano al mismo tiempo
v3 = squeeze(v(:,3,:)); %vectores v3 de los nodos
tt = [t; t; t];
v3t = (v3.*tt)';

Ks = zeros(dofsxnod*cant_puntos);
for i = 1:n
    xi = puntos(i,1);
    eta = puntos(i,2);
    zeta = puntos(i,3);

    N = [1, xi, eta, xi*eta]*A;
    Nxi = [0, 1, 0, eta]*A;
    Neta = [0, 0, 1, xi]*A;

    dN = [Nxi; Neta];

    J = [ dN*(nodos + zeta*v3t/2)
             N*(v3t)/2 ];
    dN*(nodos + zeta*v3t/2);

    invJ = J\eye(3);

    %lo pasa de (2*cant_puntos) a (3,cant_puntos)
    dN = invJ(:,1:2)*dN;

    B  = zeros(6,dofsxnod*cant_puntos);
    for inod = 1:cant_puntos
        v1 = v(:, 1, inod);
        v2 = v(:, 2, inod);

        dZN = dN(:,inod)*zeta + N(inod)*invJ(:,3);

        aux1 = [ dN(1,inod)         0          0
                         0  dN(2,inod)         0
                         0          0  dN(3,inod)
                 dN(2,inod) dN(1,inod)         0
                         0  dN(3,inod) dN(2,inod)
                 dN(3,inod)         0  dN(1,inod) ];

        aux2 = [ -v2.*dZN                        v1.*dZN
                 -v2(1)*dZN(2) - v2(2)*dZN(1)    v1(1)*dZN(2) + v1(2)*dZN(1)
                 -v2(2)*dZN(3) - v2(3)*dZN(2)    v1(2)*dZN(3) + v1(3)*dZN(2)
                 -v2(1)*dZN(3) - v2(3)*dZN(1)    v1(1)*dZN(3) + v1(3)*dZN(1) ]*0.5*t(inod);

        ini = 1 + (inod - 1)*dofsxnod;
        fin = ini + dofsxnod - 1;
        B(:,ini:fin) = [aux1 aux2];
    end%inod

    T = giro_B(J);

    B = T*B;

    mult = abs(det(J))*w(i);
    Kmin = B'*Ds*B;

    Ks = Ks + Kmin*mult;
end% i

%% Gauss Corte
[w, puntos, n] = gauss([1,1,1]);

%para hacer el jacobiano al mismo tiempo
tt = [t; t; t];
v3t = (v3.*tt)';
Kc = zeros(dofsxnod*cant_puntos);
for i = 1:n
    xi = puntos(i,1);
    eta = puntos(i,2);
    zeta = puntos(i,3);

    N = [1, xi, eta, xi*eta]*A;
    Nxi = [0, 1, 0, eta]*A;
    Neta = [0, 0, 1, xi]*A;

    dN = [Nxi; Neta];

    J = [ dN*(nodos + zeta*v3t/2)
             N*(v3t)/2 ];

    invJ = J\eye(3);

    %lo pasa de (2*cant_puntos) a (3,cant_puntos)
    dN = invJ(:,1:2)*dN; %Nx y Ny

    B  = zeros(6,dofsxnod*cant_puntos);
    for inod = 1:cant_puntos
        v1 = v(:, 1, inod);
        v2 = v(:, 2, inod);

        dZN = dN(:,inod)*zeta + N(inod)*invJ(:,3);

        aux1 = [ dN(1,inod)         0          0
                         0  dN(2,inod)         0
                         0          0  dN(3,inod)
                 dN(2,inod) dN(1,inod)         0
                         0  dN(3,inod) dN(2,inod)
                 dN(3,inod)         0  dN(1,inod) ];

        aux2 = [ -v2.*dZN                        v1.*dZN
                 -v2(1)*dZN(2) - v2(2)*dZN(1)    v1(1)*dZN(2) + v1(2)*dZN(1)
                 -v2(2)*dZN(3) - v2(3)*dZN(2)    v1(2)*dZN(3) + v1(3)*dZN(2)
                 -v2(1)*dZN(3) - v2(3)*dZN(1)    v1(1)*dZN(3) + v1(3)*dZN(1) ]*0.5*t(inod);

        ini = 1 + (inod - 1)*dofsxnod;
        fin = ini + dofsxnod - 1;
        B(:,ini:fin) = [aux1 aux2];
    end%inod

    T = giro_B(J);

    B = T*B;

    mult = abs(det(J))*w(i);
    Kmin = B'*Dc*B;

    Kc = Kc + Kmin*mult;
end% i

%% Sumar las K
K = Ks+Kc;
end



function [T] = giro_B(J)
    % J es el jacobiano
    % sistema de coordenadas local 123 en [ksi eta zeta]
    % size 5x6
    dir1 = J(1,:);
    dir3 = cross(dir1,J(2,:));
    dir2 = cross(dir3,dir1);

    M1 = [ dir1/norm(dir1); dir2/norm(dir2); dir3/norm(dir3) ];
    M2 = [ M1(:,2), M1(:,3), M1(:,1) ];
    M3 = [ M1(2,:); M1(3,:); M1(1,:) ];
    M4 = [ M2(2,:); M2(3,:); M2(1,:) ];

    T11 = M1.^2;
    T12 = M1.*M2;
    T21 = 2*M1.*M3;
    T22 = M1.*M4 + M3.*M2;

    T = [ T11 T12
          T21 T22 ];
    T(3,:) = [];
end