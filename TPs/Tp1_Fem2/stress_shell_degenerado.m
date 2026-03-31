function [Esfuerzos] = stress_shell_degenerado(nodos, Uel, E,v,t,v3)
%Stress = stress_Q4(Coord, Uel, C, Czz)
%
%Devuelve las tensiones en el punto central
%Devuelve un vector de 6 numeros
%
%Stress(1) es la tension xx
%Stress(2) es la tension yy
%Stress(3) es la tension xy
%Stress(4) es la tension xz
%Stress(5) es la tensioon zy
%
%Momento(1) es Mx
%Momento(2) es My
%Momento(3) es Mxy
%Momento(4) es Qx
%Momento(5) es Qy
%
%
%Coord es una lista (x, y) de los 4 nodos del elemento
%
%Uel son los dezplamientos de los nodos, siendo un vector
%(x1,y1,w1,alpha1,beta1,x2,y2,w2,aplha2,...) con el orden visto mas abajo
%
%E, v y t son el modulo de elasticidad, coef de poisson, espesor (1,4)
%
%poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 
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

C_b = E/(1-v^2)*[1 v 0;
                 v 1 0;
                 0 0 (1-v)/2];
C_s = c*G*[1 0; 0 1];

%% Isoparametrico
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
    aux = v3(:,i)';
    v2 = cross(aux,v1);

    v(:,:,i) = [v1;v2;aux]'; %v1 es vector fila
end%i
%% Stress en 0,0
%para hacer el jacobiano al mismo tiempo
tt = [t; t; t];
v3t = (v3.*tt)';
    
xi = 0;
eta = 0;
zeta = 1;

N = [1, xi, eta, xi*eta]*A;
Nxi = [0, 1, 0, eta]*A;
Neta = [0, 0, 1, xi]*A;

dN = [Nxi; Neta];

J = [ dN*(nodos + zeta*v3t/2)
         N*(v3t)/2 ];

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
v1 = squeeze(v(:,1,:))';
v1 = N*v1;
v2 = squeeze(v(:,2,:))';
v2 = N*v2;
v3 = N*v3';

R = [v1;v2;v3];
L(1:3,1:3) = R;
L(4:5,4:5) = R(1:2,1:2);

T = zeros(5*4);
for i=1:4
    T((i-1)*5+1:i*5, (i-1)*5+1:i*5) = L;
end
%cambiar el U a los locales, (x,y,z) da igual al de Mindlin
Uel = T*Uel;

% Esto gira B al global y lo pasa de 6x20 a 5x20 no se si se tiene que rotar al global
T2 = giro_B(J);

%B = T2*B;

    %% Tension
    z = zeta*t(1)/2; % calculamos a altura t/2 (deberiamos tmb calcularlo a -t/2)
    % OOJJJOOOOO

    eps = B*Uel
    eps = L*eps;

    S = diag([-z -z -z 1 1]);
    C = zeros(5);
    C(1:3,1:3) = C_b;
    C(4:5,4:5)= C_s;

    eps_m = eps(1:3);   % membrana
    kappa = eps(1:3);   % (si están incluidas como tales)
    gamma = eps(4:5);

    %Stress = C*S*eps; %sxx syy sxy
    Stress(1:3) = C_b*(eps_m + z*kappa);
    Stress(4:5)   = C_s*gamma;

    % Momentos(1:3) = Ds(1:3,1:3)*eps(1:3);
    % 
    % Momentos(4:5) = Dc(4:5,4:5)*eps(4:5);

    Momentos = D*eps;

    Esfuerzos(1) = Stress(1)*t(1);
    Esfuerzos(2) = Stress(2)*t(1);
    
    Esfuerzos = [Esfuerzos Momentos'];
end
%% funciones
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