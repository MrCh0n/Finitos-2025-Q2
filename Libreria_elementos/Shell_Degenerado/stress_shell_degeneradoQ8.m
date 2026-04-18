function [Esfuerzos] = stress_shell_degeneradoQ8(nodos, Uel, E,v,t,v3)
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
t_mean = mean(t);
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
cant_puntos = 8;

puntos_Q8 = [-1    -1;
              1    -1;
              1     1;
             -1     1;
              0    -1;
              1     0;
              0     1;
             -1     0];

x1 = puntos_Q8(:,1);%puntos del cuadrado
y1 = puntos_Q8(:,2);
A = [ones(cant_puntos,1) x1 y1 x1.^2 x1.*y1 y1.^2 x1.^2.*y1 y1.^2.*x1];
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
        case 5
            x1=5; x2=2; y1=5; y2=7;
        case 6
            x1=8; x2=6; y1=6; y2=3;
        case 7
            x1=7; x2=3; y1=5; y2=7;
        case 8
            x1=8; x2=6; y1=8; y2=4;
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
    
B_bot = B_local(nodos,A,v,v3t,t,0,0,-1); %bottom

B_top = B_local(nodos,A,v,v3t,t,0,0,1);

B_m = B_local(nodos,A,v,v3t,t,0,0,0);

    %% Esfuerzos
    eps_bot = B_bot*Uel; %def local

    eps_top = B_top*Uel; %eps_x, eps_y, gam_xy, gam_xz, gam_zy

    eps_m = B_m*Uel;

    Tensiones_bot =  D*eps_bot;
    Tensiones_top =  D*eps_top;
    Tensiones_m = D*eps_m;
    
    Momentos = zeros(1,7);
    Momentos(1:2) = (Tensiones_bot(1:2)+Tensiones_top(1:2))*t_mean/2;
    Momentos(3:5) = (Tensiones_bot(1:3)-Tensiones_top(1:3))*t_mean^2/12;

    dz = t_mean/2;
    [w,puntos,n] = gauss(2);
    for i = 1:n
        B = B_local(nodos,A,v,v3t,t,0,0,puntos(i));
        sigmac = D*B*Uel;
        Momentos(6:7) = Momentos(6:7)+sigmac(3:4)'*dz*w(i);
    end
    %Momentos(6:7)'-Tensiones_m(3:4)*t_mean; % era para corroborar si daba
    %lo mismo de las dos maneras
    %(Tensiones_bot(4:5)+Tensiones_top(4:5))*t_mean/2;

    Esfuerzos = Momentos;

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

function [B] = B_local(nodos,A,v,v3t,t,xi,eta,zeta)
dofsxnod = 5;
cant_puntos = 8;

N = [1, xi, eta, xi^2, xi*eta,  eta^2,  xi^2*eta,  xi*eta^2]*A;

Nxi = [0, 1, 0 2*xi, eta, 0, 2*xi*eta, eta^2]*A;%derivada de N en eta en los puntos de Gauss
Neta = [0, 0, 1, 0, xi, 2*eta, xi^2, 2*xi*eta]*A;%derivada de N en zeta

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

% Esto gira B al global y lo pasa de 6x20 a 5x20 no se si se tiene que rotar al global
T2 = giro_B(J);

B = T2*B; %esta B es en locales
end