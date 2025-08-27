function [w, puntos] = Gauss(n)
%calcula los puntos y pesos de una integral de Gauss, poner la cantidad de
%puntos en n
cant_w = ceil(n/2);
syms a [1 cant_w]
syms A [cant_w 1]


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

if flag
P(n,1) = 1;
end

%sistema a resolver
eqn = P*A == ints;

S = solve(eqn);

%% Cambiar las variables por el resultado del solve y hacer w, puntos
for i = 2:cant_w
    a(i) = S.(string(a(i)))(end);
    A(i) = S.(string(A(i)))(end);
end


A(1) = S.(string(A(1)))(end);
if ~flag
 a(1) = S.(string(a(1)))(end);
end
A = A';

%tengo los absolutos, pongo los negativos y borro el 0 extra
w = double([flip(A) A]);
puntos = double([-flip(a) a]);

if flag
 puntos(cant_w) = [];
 w(cant_w) = [];
end

end