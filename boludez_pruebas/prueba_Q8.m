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

carga = [-1 1 -1/2 1/2 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;]*-q/sqrt(2);%carga 1 es 152, 2 es 263, 3 es 374, 4 es 481
volumen = ones(8,2)*0;

R = carga_Q8(nodos,carga,volumen);

nodos = [0 0;
         1 0;
         1 1;
         0 1];

carga = [0 0 0 0;
         0 0 0 0;
         0 0 0 0;
         0 1 0 1];

volumen = [zeros(4,1) -ones(4,1)]*0;

R = carga_Q4(nodos, carga, volumen);


nodos = [0 0;
         2 0;
         1 1];

carga = [0 0 0 0;
         0 0 0 0;
         0 0 0 0];

volumen = [zeros(3,1) -ones(3,1)];

R = carga_CST(nodos, carga, volumen);

nodos = [0 0;
         1 0;
         0 1
         0.5 0;
         0.5 0.5;
         0 0.5];

carga = [0 0 0 0 0 0;
         0 10 0 1 0 1;
         0 0 0 0 0 0];

volumen = [zeros(6,1) -ones(6,1)]*0;

R = carga_LST(nodos, carga, volumen);