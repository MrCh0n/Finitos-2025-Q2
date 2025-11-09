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
