function [R] = carga_CST(nodos, carga_s, carga_v)
%Crea el vector de cargas (x,y) de un elemento CST
%
%R = carga_CST(nodos, cargas_s, cargas_v)
%
%nodos:
%   es una lista (3,2) (x, y) de los 8 nodos del elemento
%
%cargas_s:
%   son las cargas superficiales del CST, tn es la fuerza tangencial y
%   sn la normal del nodo n. Con normal saliente positiva y la tangente a su
%   derecha
%
%     ^ s
%     | ->t
%   ----- superficie 
%
%   cargas_s = [t1, s1, t2, s2;  (lado 1-2)
%               t2, s2, t3, s3;  (lado 2-3)
%               t3, s3, t1, s1;  (lado 3-1)
%
%carga_v:
%   es una lista (3,2) (x,y) de las fuerzas volumetricas en los nodos
%   ejemplo de gravedad en -y [zeros(3,1), -rho*g*t*ones(3,1)]
%   con rho la densidad del material, g la acceleracion por gravedad y t el
%   espesor
%
%poner los nodos como:
%    3
%   / \
%  /   \
% /     \
%1 - - - 2 


    %% creo el isoparametrico
    cant_puntos = 3;
    
    x1 = [0; 1; 0];
    y1 = [0; 0; 1];
    A = [ones(cant_puntos,1) x1 y1];
    A = inv(A);
    
    R = zeros(6,1);
    %% Gauss 
    %puntos de gauss normales, los transformo despues
    puntos_base = [-sqrt(3/5) 0 sqrt(3/5)];
    w_base = [5/9 8/9 5/9];
    orden = size(puntos_base,2);
    
    %puntos de las lineas alineadas con x,y
    puntos_borde = puntos_base/2+1/2;%los transformo a [0 1]
    w_borde = w_base/2;%tienen que sumar 1

    %puntos de la diagonal
    puntos_diagonal = (puntos_base+1)'/sqrt(2)*[1 -1]/sqrt(2);%los transformo a [0,sqrt(2)] y los giro -45 grados
    puntos_diagonal(:,2) = puntos_diagonal(:,2) + 1;%los muevo de {0,0}-{1 -1} a {0,1}-{1 0}
    puntos_diagonal = puntos_diagonal';%giro para comodidad
    w_diagonal = w_base/sqrt(2);%tienen que sumar sqrt(2)

    %pesos de gauss para las distintas lineas
    w = [w_borde;
        w_diagonal;
        w_borde];

    %puntos de gauss en xi para cada linea
    puntos_zeta = [zeros(orden,1)';
                    puntos_diagonal(2,:);
                    puntos_borde];
    
    %puntos de gauss en eta para cada linea 
    puntos_eta = [puntos_borde;
                   puntos_diagonal(1,:);
                   zeros(orden,1)'];
    
    nodos_linea = [1 2;
                   2 3;
                   3 1];
    
     %para sacar dx/ds, dy/ds para cada linea
     tipo = [-1 0;
              1/sqrt(2) -1/sqrt(2)
              0 1;];%direccion de tau en la linea respectiva
    
    for i = 1:3
        tau = carga_s(i,1:2:end)';
        sigma = carga_s(i,2:2:end)';
        [Rx, Ry] = superficie(nodos, A, puntos_eta(i,:), puntos_zeta(i,:), nodos_linea(i,:), tau,sigma, w(i,:) ,tipo(i,:));
        R(2*nodos_linea(i,:) - 1) = R(2*nodos_linea(i,:) - 1) + Rx;
        R(2*nodos_linea(i,:)) = R(2*nodos_linea(i,:)) + Ry;
    end
    
    %Gauss para un triangulo
    puntos =     [1/6 2/3
    1/6    1/6
    2/3    1/6];
    
    w = [1 1 1]/6;
    
    volumen =   area(nodos,puntos,w,A);
    %volumen*[1;1;1;1;1;1;1;1] la fuerza/m^2 en cada nodo, normalmente rho*g*t
    
    R(1:2:end) = R(1:2:end) + volumen*carga_v(:,1);
    R(2:2:end) = R(2:2:end) + volumen*carga_v(:,2);
end

function [F_volumen] = area(nodos,puntos,w,A)
% Gauss
orden = size(puntos,1);

F_volumen = zeros(3);

for i = 1:orden
        N = [1 puntos(i,1), puntos(i,2)]*A;%N
        Neta = [0, 1, 0]*A;%derivada de N en eta en los puntos de Gauss
        Nzeta = [0, 0, 1]*A;%derivada de N en zeta
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
        
        integrando = N'*N;
    
        F_volumen = F_volumen + integrando*det(J)*w(i);
end% i

end

function [Rx, Ry] = superficie(nodos, A,puntos_eta, puntos_zeta, nodos_linea, tau, sigma, w, direccion)
    Rx = zeros(2,1);
    Ry = zeros(2,1);
    orden = size(puntos_eta,2);
    for i = 1:orden
            N = [1 puntos_eta(i), puntos_zeta(i)]*A;%N en eta zeta =-1 
            Neta = [0, 1, 0]*A;%derivada de N en eta en los puntos de Gauss
            Nzeta = [0, 0, 1]*A;%derivada de N en zeta
            
            D = [Neta; Nzeta];
        
            J = D*nodos;
            Jeta = direccion*J;
            
            Nr = N(nodos_linea);
            Rx_aux = tau*Jeta(1)-sigma*Jeta(2);
            Ry_aux = tau*Jeta(2)+sigma*Jeta(1);
            integrando = Nr'*Nr;
    
            Rx = Rx + integrando*Rx_aux*w(i);
            Ry = Ry + integrando*Ry_aux*w(i);
    end% i
end
