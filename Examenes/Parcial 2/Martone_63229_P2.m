clc
clear
%ej 2

t = 4e-3;%m
E = 100e6;%Pa
v = 0.3;

%plain stress t<<L
C = E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];

q_max = -t*6e6;%N/m
derecha = @(x) q_max*(2000-x)/1000;
izquerda = @(x) q_max*(x)/1000;
%derecha = @(x) q_max;
%izquerda = derecha;

n=2;
x=linspace(0,2000,n+1);
for i = 1:n+1
    mesh.nodos(2*i-1,:) = [x(i) 0]; 
    mesh.nodos(2*i,:) = [x(i) 1000]; 
end

for i = 1:n
    inicial = 2*i-1;
    patron = [0 2 3 1];
    mesh.elems.con(i,:) = inicial+patron;
end

nnod = size(mesh.nodos,1);
ndof = nnod*2;
nelem = size(mesh.elems.con, 1);

K = zeros(ndof);
R = zeros(ndof,1);

mesh.dofs = 1:ndof;
mesh.dofs = reshape(mesh.dofs,2,[])';

fijos = [1 nnod-1];
movil_y = [2 nnod];

free = true(nnod, 2);

free(fijos,:) = false;
free(movil_y,1) = false;

free = reshape(free',1,[]);

for i = 1:nelem
    nodoid = mesh.elems.con(i,3:4);
    x_4 = mesh.nodos(nodoid(2),1);
    x_3 = mesh.nodos(nodoid(1),1);
    
    if x_3 < 1000
        carga3 = izquerda(x_3);
    else
         carga3 = derecha(x_3);
    end
    if x_4 < 1000
        carga4 = izquerda(x_4);
    else
         carga4 = derecha(x_4);
    end

    mesh.carga_v(i).carga = [zeros(4,1), zeros(4,1)];
    mesh.carga_s(i).carga = zeros(4,4);
    mesh.carga_s(i).carga(3,:) = [0 carga3 0 carga4];
end

% mesh.carga_s(1).carga(3,:) = [0 carga_3 0 carga_4];
% mesh.carga_s(2).carga(3,:) = [0 -q_max 0 -q_max/2];
% mesh.carga_s(2).carga(3,:) = [0 0 0 -q_max];
% mesh.carga_s(3).carga(3,:) = [0 -q_max/2 0 -q_max];

for i = 1:nelem
    nodoid = mesh.elems.con(i,:);
    coord = mesh.nodos(nodoid,:);

    Kel = t*crearK_Q4(coord, C);
    
    eledofs = reshape(mesh.dofs(nodoid,:)',1,[]);

    carga_s = mesh.carga_s(i).carga;
    carga_v = mesh.carga_v(i).carga;
    
    K(eledofs,eledofs) = K(eledofs,eledofs) + Kel;
    R(eledofs) = R(eledofs) + carga_Q4(coord, carga_s, carga_v);
end

Kr = K(free, free);
Rr = R(free);
U = zeros(ndof,1);
U(free) = Kr\Rr;

mov = reshape(U,2,[])';%cuanto se mueven

mult = 1;%cambia el multiplicador de la deformada
def = mesh.nodos+mov*mult;%deformada de los nodos

%ploteo de la deformada
plot(mesh.nodos(:,1), mesh.nodos(:,2),'*')
hold on
for i = 1:nelem
    nodoid = mesh.elems.con(i,:);
    grafico = [def(nodoid,:);def(nodoid(1),:)];
    plot(grafico(:,1), grafico(:,2))
end

mesh.elems.stress = zeros(nnod,5);%tension xx yy xy y principales el indice es el nodo al que pertenece
mesh.elems.stress_gauss = zeros(4*nelem,3);
cont = zeros(nnod,1);
for i = 1:nelem
    nodoid = mesh.elems.con(i,:);
    gauss_id = [1:4] + 4*(i-1);

    nodos = mesh.nodos(nodoid,:);

    dir = mesh.dofs(nodoid,:);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas

    Uel = U(dir);
    [Stress_nodos, Stress_Gauss] = stress_Q4(nodos,Uel,C);
    mesh.elems.stress(nodoid,:) = mesh.elems.stress(nodoid,:)+Stress_nodos;
    mesh.elems.stress_gauss(gauss_id,:) = Stress_Gauss;
    cont(nodoid) = cont(nodoid) +1;
end
mesh.elems.stress = mesh.elems.stress./cont;
%Imprimir resultados
fprintf('\n          \tTensiones nodales (MPa)      \n')
fprintf('\t Sigma xx\t Sigma yy\t Sigma xy\tSigma I\t Sigma II\n')
for i=1:nnod
    fprintf('nodo %d\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n', i, mesh.elems.stress(i,:)/1e6);
end

fprintf('\n          \tTensiones nodales Gauss (MPa)      \n')
fprintf('\t\t\t Sigma xx\t Sigma yy\t Sigma xy\n')
for i=1:1
    fprintf('\nelemento %d\n',i);
    offset = 4*(i-1);
    for j= 1:4
        fprintf('\tnodo gauss %d\t %.5f\t %.5f\t %.5f\n', j, mesh.elems.stress_gauss(offset+j,:)/1e6);
    end
end

function [R] = carga_Q4(nodos, carga_s, carga_v, options)
%Crea el vector de cargas (x,y) de un elemento Q4
%
%R = carga_Q4(nodos, cargas_s, cargas_v)
%
%nodos:
%   es una lista (4,2) (x, y) de los 8 nodos del elemento
%
%cargas_s:
%   son las cargas superficiales del Q4, tn es la fuerza tangencial y
%   sn la normal del nodo n. Con normal saliente positiva y la tangente a su
%   derecha
%
%     ^ s
%     | ->t
%   ----- superficie 
%
%   cargas_s = [t1, s1, t2, s2;  (lado 1-2)
%               t2, s2, t3, s3;  (lado 2-3)
%               t3, s3, t4, s4;  (lado 3-4)
%               t4, s4, t1, s1]  (lado 4-1)
%
%carga_v:
%   es una lista (4,2) (x,y) de las fuerzas volumetricas en los nodos
%   ejemplo de gravedad en -y [zeros(4,1), -rho*g*t*ones(4,1)]
%   con rho la densidad del material, g la acceleracion por gravedad y t el
%   espesor
%
% poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 

arguments
        nodos (4,2) {mustBeNumeric}
        carga_s (4,4) {mustBeNumeric}
        carga_v (4,2) {mustBeNumeric}
        options.coordenada (1,1) {mustBeNumericOrLogical} = false
    end



    %% creo el isoparametrico
    cant_puntos = 4;
    
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(cant_puntos,1) x1 y1 x1.*y1];
    A = inv(A);
    
    R = zeros(8,1);

    %para sacar dx/ds, dy/ds para cada linea
     tipo = [-1 0;
              0 -1;
              1 0;
              0 1;];%direccion de tau en la linea respectiva

     %giro la carga para que sea en x y
     if options.coordenada

        puntos_zeta = [y1(1:2)';
                        y1(2:3)';
                        y1([3 4])';
                        y1([4 1])'];
        
        puntos_eta = [x1(1:2)';
                        x1(2:3)';
                        x1(3:4)';
                        x1([4 1])'];

        posicion = [1 2;
                    3 4];

        for i = 1:4
            carga_s(i,:) = girar_carga(nodos, carga_s(i,:), A, puntos_eta(i,:), puntos_zeta(i,:), tipo(i,:), posicion);
        end
     end
    %% Gauss
    puntos = [-sqrt(3/5) 0 sqrt(3/5)];
    w = [5/9 8/9 5/9];
    
    orden = size(puntos,2);
    
    puntos_zeta = [-ones(orden,1)';
                    puntos;
                    ones(orden,1)';
                    puntos];
    
    puntos_eta = [puntos;
                   ones(orden,1)';
                   puntos
                   -ones(orden,1)'];
    
    nodos_linea = [1 2;
                   2 3;
                   3 4;
                   4 1];
    
    for i = 1:4
        tau = carga_s(i,1:2:end)';
        sigma = carga_s(i,2:2:end)';
        [Rx, Ry] = superficie(nodos, A, puntos_eta(i,:), puntos_zeta(i,:), nodos_linea(i,:), tau,sigma,w,tipo(i,:));
        R(2*nodos_linea(i,:) - 1) = R(2*nodos_linea(i,:) - 1) + Rx;
        R(2*nodos_linea(i,:)) = R(2*nodos_linea(i,:)) + Ry;
    end
    
    volumen =   area(nodos,puntos,w,A);
    %volumen*[1;1;1;1;1;1;1;1] la fuerza/m^2 en cada nodo, normalmente rho*g*t
    
    R(1:2:end) = R(1:2:end) + volumen*carga_v(:,1);
    R(2:2:end) = R(2:2:end) + volumen*carga_v(:,2);
end

function [F_volumen] = area(nodos,puntos,w,A)
% Gauss
orden = size(puntos,2);

F_volumen = zeros(4);

for i = 1:orden
    for j = 1:orden
        N = [1 puntos(i), puntos(j), puntos(i)*puntos(j)]*A;%N
        Neta = [0, 1, 0, puntos(j)]*A;%derivada de N en eta en los puntos de Gauss
        Nzeta = [0, 0, 1, puntos(i)]*A;%derivada de N en zeta
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
        
        integrando = N'*N;
    
        F_volumen = F_volumen + integrando*det(J)*w(i)*w(j);
    end% j
end% i

end

function [Rx, Ry] = superficie(nodos, A,puntos_eta, puntos_zeta, nodos_linea, tau, sigma, w, direccion)
    Rx = zeros(2,1);
    Ry = zeros(2,1);
    orden = size(puntos_eta,2);
    for i = 1:orden
            N = [1 puntos_eta(i), puntos_zeta(i), puntos_eta(i)*puntos_zeta(i)]*A;%N en eta zeta =-1 
            Neta = [0, 1, 0,  puntos_zeta(i)]*A;%derivada de N en eta en los puntos de Gauss
            Nzeta = [0, 0, 1, puntos_eta(i)]*A;%derivada de N en zeta
            
            D = [Neta; Nzeta];
        
            J = D*nodos;
            d_x_y = direccion*J;%[dx/ds dy/ds]
            
            Nr = N(nodos_linea);
            Rx_aux = tau*d_x_y(1)-sigma*d_x_y(2);
            Ry_aux = tau*d_x_y(2)+sigma*d_x_y(1);
            integrando = Nr'*Nr;
    
            Rx = Rx + integrando*Rx_aux*w(i);
            Ry = Ry + integrando*Ry_aux*w(i);
    end% i
end


function [carga_s] = girar_carga(nodos, carga_s, A, puntos_eta, puntos_zeta, direccion,nodo_linea)
    orden = size(puntos_eta,2);
    for i = 1:orden
            Neta = [0, 1, 0,  puntos_zeta(i)]*A;%derivada de N en eta en los puntos de Gauss
            Nzeta = [0, 0, 1, puntos_eta(i)]*A;%derivada de N en zeta
            
            D = [Neta; Nzeta];
        
            J = D*nodos;
            tangencial = direccion*J;%[dx/ds dy/ds]
            normal = tangencial*[0 1;-1 0];
            giro = [tangencial' normal']/norm(tangencial);
            
            carga_s(nodo_linea(i,:)) =  carga_s(nodo_linea(i,:))*giro;
            carga_s(nodo_linea(i,:))*giro;
    end% i
end


function [K] = crearK_Q4(nodos,C)
%Crea la matriz de resistencia de un Q4 sin espesor
%
%K = crearK_Q4(nodos, C)
%
%nodos es una lista (x, y) de los 4 nodos del elemento
%C es el constitutivo de 3x3
%
%poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 

arguments (Input)
   nodos(4,2) {mustBeNumeric}
   C(3,3) {mustBeNumeric}
end
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
    
        J = D*nodos;
    
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
        
        K = K + Kmin*mult;

    end% j
end% i

end

function [Stress, sigma_rs] = stress_Q4(Coord, Uel, C)
    Stress = zeros(4,5); %[xx yy xy principales]
    
    % isoparametrico
    cant_puntos = 4;
    
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(cant_puntos,1) x1 y1 x1.*y1];
    
    %iso para extrapolar las tensiones en puntos de Gauss/superonvergentes
    x2 = [-1; -1; 1; 1];
    y2 = [-1; 1; -1; 1];
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