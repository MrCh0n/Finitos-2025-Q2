function [] = plot_CST(nodos, elem, dofs, U, mult, C)

nelem = size(elem,1);

plot(nodos(:,1), nodos(:,2), 'k*');
hold on

todo = mult*reshape(U,2,[])';
deformada = nodos + todo;
plot(deformada(:,1), deformada(:,2), 'bo');

for i = 1:nelem
    dir_nod = elem(i, :);

    puntos = nodos(dir_nod,:);

    x = puntos(:,1);
    y = puntos(:,2);
    one = ones(size(puntos,1),1);
    A = [one x y];
    
    Bx = [0 1 0]/A;
    By = [0 0 1]/A;
    
    B(1, [1 3 5]) = Bx;
    B(2, [2 4 6]) = By;
    B(3, :) = [By(1) Bx(1) By(2) Bx(2) By(3) Bx(3)];
    
    dir = reshape(dofs(dir_nod,:)', [], 1);

    Uloc = U(dir);
    Stress(i,:) = C*B*Uloc;
end

ncol = 1024;
cmap = jet(ncol);
colormap(cmap);
for i =1:3
    Stressmax(i) = max(Stress(:,i));
    Stressmin(i) = min(Stress(:,i));
end

Ylabel = ['Sress XX [MPa]';
        'Sress YY [MPa]';
        'Sress XY [MPa]'];
hold off
for j = 1
    figure(j);

    for i = 1:nelem
        nodo1 = elem(i, 1);
        nodo2 = elem(i, 2);
        nodo3 = elem(i, 3);
        
        triangulo = [deformada(nodo1, :); deformada(nodo2, :); deformada(nodo3, :)];
        
        T = (Stress(i,j)-Stressmin(j))/(Stressmax(j)-Stressmin(j));
        idx = round(1 + T*(ncol-1));
        color = cmap(idx,:);
    
        plot(polyshape(triangulo), 'FaceColor', color);
        hold on
    end
clim([Stressmin(j) Stressmax(j)]/1e6);
ylabel(colorbar, Ylabel(j,:))
hold off
end

end