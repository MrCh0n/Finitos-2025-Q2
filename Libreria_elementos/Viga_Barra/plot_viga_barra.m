function [] = plot_viga_barra(nodos,elem,dofs, U, divisiones, mult)

nelem = size(elem,1);

Nloc = zeros(2*divisiones,6);
plot(nodos(:,1), nodos(:,2), 'k*');
hold on

todo = mult*reshape(U,3,[])';
deformada = nodos + todo(:, [1 2]);
plot(deformada(:,1), deformada(:,2), 'bo');
for i = 1:nelem
    nodo1 = elem(i,1);
    nodo2 = elem(i,2);

    V = nodos(nodo2,:) - nodos(nodo1,:);%vector de direccion
    L = norm(V);

    %cosenos directores
    cs = V/L;
    c = cs(1);
    s = cs(2);
  
    %matriz de rotacion
    Q = [c s 0;
         -s c 0;
         0 0 1];
    T(1:3, 1:3) = Q;
    T(4:6, 4:6) = Q;

    dir = [dofs(nodo1,:) dofs(nodo2,:)];
   
    Uloc =  T*U(dir);
    
    %creo la matrix de N    
    x = linspace(0,L, divisiones);
    
    dir_bar = [1 4];
    dir_vig = [2 3 5 6];
    
    Nel_bar = [L-x; x]'/L;
    Nel_vig = [(L^3 - 3*L*x.^2 + 2*x.^3)/L^3;
               (L^2*x - 2*L*x.^2 + x.^3)/L^2;
               (- 2*x.^3 + 3*L*x.^2)/L^3;
               -(- x.^3 + L*x.^2)/L^2]';

    Nloc(1:2:2*divisiones,dir_bar) = Nel_bar;
    Nloc(2:2:2*divisiones,dir_vig) = Nel_vig;
    
    %crea los desplazamientos u,v locales del elemento
    %mult amplifica la deformacion
    Udiv = reshape(Nloc*Uloc,2,divisiones)'*mult;
    
    %para crear los puntos iniciales de las divisiones
    inicial = nodos(nodo1,:);
    final = nodos(nodo2,:);
    deltaL = (final-inicial)/(divisiones-1);
    
    Q = [c s;
         -s c];%para girar al global
    for j = 1:divisiones
        Udiv(j,:) = (Q'*Udiv(j,:)')' + inicial + deltaL*(j-1); %lo giro al global y le sumo su punto inicial
    end

    plot(Udiv(:,1), Udiv(:,2))
end
hold off

end