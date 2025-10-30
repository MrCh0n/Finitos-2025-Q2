clc
clear

E = 200e9;%modulo de young
v = 0.3;%poisson
alpha = 12e-6;%dilatacion termica

Pi = 25e6;
interferencia_radial = 0.5e-3;%m

C = E/((1+v)*(1-2*v))*[1-v v 0;v 1-v 0;0 0 (1-2*v)/2];%sin t

ri = 250e-3;%radio interno m
re = 500e-3;%radio externo m

dt = -interferencia_radial/re/alpha;%seria la temperatura equivalente

grados_45 = [1 1]/sqrt(2);

mesh.nodos = [ri 0;
             ri*grados_45;
             0 ri;
             0 re;
             re*grados_45;
             re 0];

mesh.elems.con = [1 6 5 2;
                  2 5 4 3];

nnod = size(mesh.nodos,1);
ndof = 2*nnod;%tienen 2 dof
nelem = size(mesh.elems.con,1);

mesh.dofs = reshape(1:ndof,2,[])';

mesh.free = true(nnod,2);
movil_y = [3 4];
mesh.free(movil_y,1) = false;
movil_x = [1 6];
mesh.free(movil_x,2) = false;

mesh.free = reshape(mesh.free',[],1);

R = zeros(ndof,1);

Pi = Pi*pi/4*ri;

R(mesh.dofs(1,:)) = R(mesh.dofs(1,:)) + Pi*[1; 0];
R(mesh.dofs(2,:)) = R(mesh.dofs(2,:)) + Pi*grados_45';
R(mesh.dofs(3,:)) = R(mesh.dofs(3,:)) + Pi*[0; 1];
% 
% for i = 1:nelem%esto es como lo haria si funcionara cargas Q4
%     nodoid = mesh.elems.con(i,:);
% 
%     nodos = mesh.nodos(nodoid,:);
% 
%     carga = cargas_Q4(nodos, [1 4]);
% 
%     dir = mesh.dofs(nodoid,:);
%     dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas
% 
%     Fx = carga*R([dir(1) dir(7)]);
%     Fy = carga*R([dir(2) dir(8)]);
%     R(dir) = [Fx(1); Fy(1); Fx(2); Fy(2)];
% end

K = zeros(ndof);


for i = 1:nelem
    nodoid = mesh.elems.con(i,:);

    nodos = mesh.nodos(nodoid,:);

    Kel = crearK_Q4(nodos,C);

    dir = mesh.dofs(nodoid,:);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas
        
    K(dir,dir) = K(dir,dir) + Kel;
end

Ft = alpha*dt*E*pi/4*re;

R(mesh.dofs(4,:)) = R(mesh.dofs(4,:)) + Ft*[0; -1];
R(mesh.dofs(5,:)) = R(mesh.dofs(5,:)) + Ft*grados_45'*-1;
R(mesh.dofs(6,:)) = R(mesh.dofs(6,:)) + Ft*[-1; 0];%presion termica en los nodos, por lo mismo que la otra R

Kr = K(mesh.free,mesh.free);

Rr = R(mesh.free);

U = zeros(ndof,1);
U(mesh.free) = Kr\Rr;

stress = zeros(nnod,5);%tension xx yy xy y principales el indice es el nodo al que pertenece
for i = 1:nelem
    nodoid = mesh.elems.con(i,:);

    nodos = mesh.nodos(nodoid,:);

    dir = mesh.dofs(nodoid,:);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas

    Uel = U(dir);
        
    stress(nodoid,:) = stress_Q4(nodos,Uel,C);
end


function [K] = crearK_Q4(nodos,C)
%% creo el isoparametrico
cant_puntos = 4;

x1 = [-1; 1; 1; -1];
y1 = [-1; -1; 1; 1];
A = [ones(cant_puntos,1) x1 y1 x1.*y1];


%% Gauss
puntos = [-sqrt(3/5) 0 sqrt(3/5)];
w = [5/9 8/9 5/9];

orden = size(puntos,2);

B = zeros(3,2*cant_puntos);
K = 0;

for i = 1:orden
    for j = 1:orden
        Neta = [0, 1, 0 puntos(j)]/A;
        Nzeta = [0, 0, 1, puntos(i)]/A;
        
        D = [Neta; Nzeta];
    
        J = D*nodos;%jacobiano
    
        Bs = J\D;
            
        Bx = Bs(1,:);
        By = Bs(2,:);
        
        dir1 = 1:2:2*cant_puntos;
        dir2 = 2:2:2*cant_puntos;
        %crear la matrz B
        B(1,dir1) = Bx;
        B(2,dir2) = By;
        B(3,dir1) = By;
        B(3,dir2) = Bx;
    
        mult = abs(det(J))*w(i)*w(j);
        Kmin = B'*C*B;
        
        K = K + Kmin*mult;%integral
    end% j
end% i

end

function [cargas] = cargas_Q4(nodos, elejidos)
    %% creo el isoparametrico
    cant_puntos = 4;
    
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(cant_puntos,1) x1 y1 x1.*y1];

    syms eta xi
    N = [1 eta xi xi*eta]/A;

    Nreducido = N(elejidos);

    dx = diff(Nreducido,eta)*nodos(elejidos,1);
    dy = diff(Nreducido,xi)*nodos(elejidos,2);

    ds = sqrt(dx^2+dy^2);

    integrando = Nreducido*ds;

    if x1(elejidos(1)) == x1(elejidos(2)) 
        integrando = subs(integrando,eta,x1(elejidos(1)));
        cargas = int(integrando,-1,1);
    else
        integrando = subs(integrando,xi,y1(elejidos(1)));
        cargas = int(integrando,-1,1);
    end
end


function [Stress] = stress_Q4(Coord, Uel, C)
    Stress = zeros(4,5); %[xx yy xy principales]
    
    % isoparametrico
    cant_puntos = 4;
    
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(cant_puntos,1) x1 y1 x1.*y1];
    
    %iso para extrapolar las tensiones en puntos de Gauss/superonvergentes
    x2 = [-1; 1; 1; -1];
    y2 = [-1; -1; 1; 1];
    A_rs = [ones(4,1) x2 y2 x2.*y2];
    
    r = x1*sqrt(3);%las esquinas de xi eta
    s = y1*sqrt(3);

    %donde quiero la tensiones
    puntos = [-sqrt(1/3) sqrt(1/3)];
     %toda la seccion superior se cambia junta

    dir1 = 1:2:2*cant_puntos;
    dir2 = 2:2:2*cant_puntos;
    
    cant = size(puntos,2);
    sigma_rs = zeros(cant^2,3);% (:,1) xx, 2 yy, 3 xy, 4 zz
    
    %tension en los puntos de Gauss
    for i = 1:cant
        offset = (i-1)*cant;
        for j = 1:cant
            Neta = [0, 1, 0 puntos(j)]/A;
            Nzeta = [0, 0, 1, puntos(i)]/A;
                
            D = [Neta; Nzeta];
        
            J = D*Coord;
        
            Bs = J\D;
                
            Bx = Bs(1,:);
            By = Bs(2,:);
        
            %crear la matrz B
            Bel(1,dir1) = Bx;
            Bel(2,dir2) = By;
            Bel(3,dir1) = By;
            Bel(3,dir2) = Bx;
        
            sigma_rs(offset+j,1:3) = C*Bel*Uel; %sxx syy sxy
        end
    end
    
    for j = 1:cant_puntos
        N = [1 r(j) s(j) r(j)*s(j)]/A_rs;
        sigma = N*sigma_rs;%las tensiones en los nodos del elemento
 

        Stress(j,1:3) = sigma;
    
        sxx = Stress(j,1);
        syy = Stress(j,2);
        sxy = Stress(j,3);
    
        sigma_plano = [sxx sxy;
                       sxy syy];
    
        Stress(j,4:5) = eig(sigma_plano);
    end
end