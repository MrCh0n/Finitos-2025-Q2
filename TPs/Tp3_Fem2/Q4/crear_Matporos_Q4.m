function [K,Cg,F,Mm,Kp] = crear_Matporos_Q4(nodos,elems,dofs,C,t,alpha,M,k,mu_f,dt)
dofselem = 8;
nnod = size(nodos,1);
nelem = size(elems,1);
ndof = numel(nodos);
nnz = nelem*dofselem^2; %si ningun nodo se repite se tienen esta cantidad de posibles no zeros

I = zeros(nnz,1);
J = zeros(nnz,1);
V = zeros(nnz,1);
cont = 1;

Cg = zeros(ndof, nnod);
F = zeros(nnod, ndof);
Mm = zeros(nnod, nnod);
Kp = zeros(nnod, nnod);
for i = 1:nelem
    nodoid = elems(i,:);

    coord = nodos(nodoid,:);

    Kel = t*crearK_Q4(coord, C);

    Cg_el = crearCg_Q4(coord, t, alpha);

    F_el = crearF_Q4(coord, t, alpha, dt);

    M_el = crearM_Q4(coord, M, t);

    Kp_el = crearK_Q4_poros(coord, k, mu_f, t);

    dir = reshape(dofs(nodoid,:)',1,[]);

    Cg(dir,nodoid) = Cg(dir,nodoid) + Cg_el;
    F(nodoid,dir) = F(nodoid,dir) + F_el;
    Mm(nodoid, nodoid) = Mm(nodoid, nodoid) + M_el;
    Kp(nodoid, nodoid) = Kp(nodoid, nodoid) + Kp_el;

    for a = 1:dofselem
        for b = 1:dofselem
            I(cont) = dir(a);
            J(cont) = dir(b);
            V(cont) = Kel(a,b);
            cont = cont+1;
        end%b
    end%a
end
K = sparse(I,J,V);
end