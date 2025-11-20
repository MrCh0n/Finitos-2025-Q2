clc
clear

a = 1;
q = 1;

nodos = [2*a 0;
        0 2*a;
        0 a;
        a 0;
        a a;
        0 a*3/2
        a/2 a/2
        a*3/2 0];

carga = [0 1 0 1/2 0 0;
        0 0 0 0 0 0;
        1 0 1 0 1 0;
        0 0 0 0 0 0;]*-q;%carga 1 es 152, 2 es 263, 3 es 374, 4 es 481
volumen = ones(8,2)*0;

R = carga_Q8(nodos,carga,volumen,'coordenada', true);
R/2/sqrt(2);

nodos = [0 0;
         1 0;
         0.5 1;
         0 1];

carga = [0 0 0 0;
         0 0  0 0;
         0 5860 0 5860;
         0 0 0 0]*0;

volumen = [zeros(4,1) -2400*9.81*ones(4,1)];
volumen = [zeros(4,1) -ones(4,1)];

R = carga_Q4(nodos, carga, volumen, 'coordenada', true)


E = 200e9;
t = 1;
v = 0.3;


C = t*E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];
K = crearK_Q4(nodos, C);


libres = 1:8;
libres([1 2 3 4]) = [];
Kr = K(libres, libres);
Rr = R(libres);
U = zeros(8,1);
U(libres) = Kr\Rr;

mov = reshape(U,2,[])';

def = nodos+mov/max(abs(U))/2e1;

plot(def(:,1), def(:,2))
hold on
plot(nodos(:,1), nodos(:,2))


nodos = [0 0;
         1 0;
         0 1];

carga = [0 0 0 0;
         0 1 0 1;
         0 0 0 0];
 
volumen = [zeros(3,1) -ones(3,1)]*0;

R = carga_CST(nodos, carga, volumen, 'coordenada', true);
R = carga_CST(nodos, carga, volumen);

nodos = [0 0;
         1 0;
         0 1
         0.5 0;
         0.5 0.5;
         0 0.5];

carga = [0 0 0 0 0 0;
         0 1 0 1 0 1;
         0 0 0 0 0 0];

volumen = [zeros(6,1) -ones(6,1)]*0;

R = carga_LST(nodos, carga, volumen);
R = carga_LST(nodos, carga, volumen, 'coordenada',true);
