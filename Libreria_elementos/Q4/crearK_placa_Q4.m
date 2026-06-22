function [K] = crearK_placa_Q4(nodos, E, v, t)
%Crea la matriz de resistencia de una placa Q4
%
%K = crearK_placa_Q4(nodos, E, v, t)
%
%nodos es una lista (x, y) de los 4 nodos del elemento
%
%E es el modulo de young del material
%v el modulo de poisson del material
%t el espesor de la placa
%
%poner los nodos como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 
%% placas de Kirchhoff
    C = E*t^3/(12*(1-v^2))*[1 v 0; v 1 0; 0 0 (1-v)/2];

%% Creo el A del isoparametrico placa y para el jacobiano
cant_nodos = 4;

x = [-1; 1; 1; -1];
y = [-1; -1; 1; 1];

%para el cambio de variable
A_jac = [ones(cant_nodos,1) x y x.*y];
A_jac = inv(A_jac);


% 3 dofs/nodo
% 12 funcionalidades en total, 4 para el normal, 4 para d/dx y 4 para d/dy
% w : [1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x*y^3]
% tita_z(x) : [0 1 0 2*x y 0 3x^2 2*x*y y^2 0 3*x^2*y y^3]
% tita_z(y) : [0 0 1 0 x 2*y 0 x^2 2*x*y 3*y^2 x^3 3*x*y^2]
cant_puntos = cant_nodos*3;

uno = ones(cant_nodos,1);
zero = zeros(cant_nodos,1);

A_w = [uno x y x.^2 x.*y y.^2 x.^3 x.^2.*y x.*y.^2 y.^3 x.^3.*y x.*y.^3];
A_tita_x = [zero uno zero 2*x y zero 3*x.^2 2*x.*y y.^2 zero 3*x.^2.*y y.^3];
A_tita_y = [zero zero uno zero x 2*y zero x.^2 2*x.*y 3*y.^2 x.^3 3*x.*y.^2];

%para las funciones de forma
A = zeros(cant_puntos);
for i = 1:cant_nodos
    aux = (1:3) + (i-1)*3;

    A(aux,:) = [A_w(i,:); A_tita_x(i,:); A_tita_y(i,:)];
end
A = inv(A); % para no invertirla siempre en el for loop de Gauss
% syms x y
% 
% N = [1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x*y^3]*A;
% diff(N(1),x)
% reshape(N,[],1)
% subs(N,[x y], [0 0])

%% Gauss
puntos = [-sqrt(3/5) 0 sqrt(3/5)];
w = [5/9 8/9 5/9];

orden = size(puntos,2);

B = zeros(3,cant_puntos);
K = 0;

for i = 1:orden
    for j = 1:orden
        %derivadas parciales en el isoparemetrico
        Nxx = [0 0 0 2 0 0 6*puntos(i) 2*puntos(j) 0 0 6*puntos(j)*puntos(i) 0]*A;
        Nyy = [0 0 0 0 0 2 0 0 2*puntos(i) 6*puntos(j) 0 6*puntos(i)*puntos(j)]*A;
        Nxy = [0 0 0 0 1 0 0 2*puntos(i) 2*puntos(j) 0 3*puntos(i)^2 3*puntos(j)^2]*A;
        

        Nxxx = [0 0 0 0 0 0 6 0 0 0 6*puntos(j) 0]*A;
        Nyyy = [0 0 0 0 0 0 0 0 0 6 0 6*puntos(i)]*A;

        Nxxy = [0 0 0 0 0 0 0 2 0 0 6*puntos(i) 0]*A; 
        Nyyx = [0 0 0 0 0 0 0 0 2 0 0 6*puntos(j)]*A;

        
        Neta = [0, 1, 0 puntos(j)]*A_jac;
        Nzeta = [0, 0, 1, puntos(i)]*A_jac;
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
    
        % Bs = J\D;
        % 
        % Bx = Bs(1,:);
        % By = Bs(2,:); 
        B_w = ones(3,4);
        B_tita_x = 2*ones(3,4);
        B_tita_y = 3*ones(3,4);

        dir_w = 1:3:cant_puntos;
        dir_x = 2:3:cant_puntos;
        dir_y = 3:3:cant_puntos;

        % %crear la matrz B
        B(:,dir_w) = B_w;
        B(:,dir_x) = B_tita_x;
        B(:,dir_y) = B_tita_y;

        %frr = fxx * cos^2 + fxy sen*cos + fyy sen^2
  
        mult = abs(det(J))*w(i)*w(j);
        Kmin = B'*C*B;
        
        K = K + Kmin*mult;

    end% j
end% i


end