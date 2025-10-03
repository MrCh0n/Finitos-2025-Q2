clc
clear

syms x y x1 x2 x3 y1 y2 y3

%% Triangulo cte
nodo1 = [x1 y1];
nodo2 = [x2 y2];
nodo3 = [x3 y3];

X = [1; x; y];

A = [1 nodo1;
    1 nodo2;
    1 nodo3];

N = A\X;

triangulo = [nodo1;nodo2;nodo3];

%plot(polyshape(triangulo), 'Facecolor', 'g')

%% Cuadrado
%Q4
nodo1 = [0.5 0.5];
nodo2 = [-0.5 0.5];
nodo3 = [0.5 -0.5];
nodo4 = [-0.5 -0.5];

X = [1 x y x*y];

A = [subs(X, [x y], nodo1);
    1 nodo2 nodo2(1)*nodo2(2);
    1 nodo3 nodo3(1)*nodo3(2);
    1 nodo4 nodo4(1)*nodo4(2)];

N = X/A
cant_puntos = 4;

Bx = diff(N,x);
By = diff(N,y);

        dir1 = 1:2:2*cant_puntos;
        dir2 = 2:2:2*cant_puntos;
        %crear la matrz B
        B(1,dir1) = Bx;
        B(2,dir2) = By;
        B(3,dir1) = By;
        B(3,dir2) = Bx;
%% Triangulo con Q4
x1 = [-1;1;-1;1];
y1 = [-1;-1;1;1];
A = [ones(4,1) x1 y1 x1.*y1];

N = [1 x y x*y]/A;

triangulo = [0 0; 1 0;0 1];

mapeo = [triangulo(1,:); triangulo(2,:); triangulo(3,:); triangulo(3,:)];
N2 = N*mapeo;

puntos = [-1/sqrt(3) 1/sqrt(3)];

J = [diff(N,x); diff(N,y)]*mapeo;

u = N2(1);
v = N2(2);

integrando = u^3*det(J);

Area = 0;
for i = 1:2
    for j = 1:2
        idx = 2*(i-1)+j;
        nodos(idx,:) = subs(N2, [x y], [puntos(i), puntos(j)]);
        Area = Area + subs(integrando, [x y], [puntos(i) puntos(j)]);
    end
end


% plot(nodos(:,1), nodos(:,2),'xk')
% hold on
% 
plot(polyshape(triangulo))
hold on
% hold off

%% triangulo FST bien
x1 = [0; 1; 0; 0.5; 0; 0.5];
y1 = [0; 0; 1; 0; 0.5; 0.5];
A = [ones(6,1) x1 y1 x1.*y1 x1.^2 y1.^2];

N = ([1 x y y*x x^2 y^2]/A);

puntos =     [1/6 2/3
    1/6    1/6
    2/3    1/6];

w = [1 1 1]/6;

triangulo = [0 0; 2 0; 0 2; 0.5 0; 0 0.5; 0.5 0.5];

N2 = N*triangulo;

J = [diff(N,x); diff(N,y)]*triangulo;

u = N2(1);
v = N2(2);

integrando = det(J);

Area = 0;
for i = 1:3
        nodos(i,:) = subs(N2, [x y], [puntos(i,:)]);
        Area = Area + w(i)*subs(integrando, [x y], [puntos(i,:)]);
end


plot(nodos(:,1), nodos(:,2),'xk')
hold on

% plot(polyshape(triangulo))
% hold off
