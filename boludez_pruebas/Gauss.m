function [w, puntos] = Gauss(n)
%calcula los puntos y pesos de una integral de Gauss, poner la cantidad de
%puntos en n
cant_w = ceil(n/2);
syms a [1 cant_w]
syms A [cant_w 1]
b = a;

flag = mod(n,2) == 1;%si hay un flag es porque el 0 es un caso especial

if flag
a(1) = 0;
end

ints = zeros(n,1);

for i = 1:n-1
    P(i,:) = a.^(2*i);%fila de potencias
    ints(i) = 1/(2*i+1);%integrales de pares -1 1
end

ints(n) = 2;%la suma de pesos tienen que dar 2
P(n,:) = 2*ones(cant_w,1);

idx1 = 1:cant_w;
idx2 = cant_w+1:n;

if flag
P(n,1) = 1;
idx1 = 2:cant_w;
idx2 = [1 cant_w+1:n];
end



A(idx1) = P(idx1,idx1)\ints(idx1);
%sistema a resolver
eqn = P(idx2,1:cant_w)*A == ints(idx2);

S = solve(eqn);
%% Cambiar las variables por el resultado del solve y hacer w, puntos
for i = 2:cant_w
    a(i) = abs(S.(string(a(i)))(end));
end

if ~flag
 a(1) = abs(S.(string(a(1)))(end));
else
    A(1) = abs(S.(string(A(1)))(end));
end

A = subs(A,b,a);
A = A';

%tengo los absolutos, pongo los negativos y borro el 0 extra
w = double([flip(A) A]);
puntos = double([-flip(a) a]);

if flag
 puntos(cant_w) = [];
 w(cant_w) = [];
end

end
