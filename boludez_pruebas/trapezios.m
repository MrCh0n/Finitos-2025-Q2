% clc
% clear

nodos = [0 0;
         1 0;
         1 1;
          0 1;
          2 0;
          2 1];

elems = [1 2 3 4;
         2 5 6 3];

nnod = size(nodos, 1);
ndof = 2*nnod;
nelem = size(elems, 1);

carga_s = [0 0 0 0;
           0 0 0 0;
           0 -1 0 -1;
           0 0 0 0];

carga_v = [zeros(4,1), -ones(4,1)]*0;

K = zeros(ndof);
R = zeros(ndof,1);

dofs = 1:ndof;
dofs = reshape(dofs,2,[])';

empotrados = [1 2 5];
free = true(nnod, 2);

free(empotrados,:) = false;

free = reshape(free',1,[]);

E = 200e9;
t = 1;
v = 0.3;


C = t*E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];

for i = 1:nelem
    nodoid = elems(i,:);
    coord = nodos(nodoid,:);

    Kel = crearK_Q4(coord, C);
    
    eledofs = reshape(dofs(nodoid,:)',1,[]);
    
    K(eledofs,eledofs) = K(eledofs,eledofs) + Kel;
    R(eledofs) = R(eledofs) + carga_Q4(coord, carga_s, carga_v,"coordenada",true);
end

Kr = K(free, free);
Rr = R(free);
U = zeros(8,1);
U(free) = Kr\Rr;

mov = reshape(U,2,[])';

def = nodos+mov*5e10;

plot(nodos(:,1), nodos(:,2),'*')
hold on
grafico = [def(elems(1,:),:);def(1,:);def(elems(2,:),:)];
plot(grafico(:,1), grafico(:,2))