function [val] = check_Q4(nodos)
    distorcion = distortion(nodos);
    stretch = alarge(nodos);

    dis = distorcion>0.6;
    stretchs = stretch>0.3;

    val = dis && stretchs;
end

function [dis] = distortion(nodos)
%% creo el isoparametrico
cant_puntos = 4;

x1 = [-1; 1; 1; -1];
y1 = [-1; -1; 1; 1];
A = [ones(cant_puntos,1) x1 y1 x1.*y1];
A=inv(A);


%% Gauss
[w, puntos, n] = gauss([2,2]);

Area = 0;

for i = 1:n
        xi = puntos(i,1);
        eta = puntos(i,2);

        Neta = [0, 1, 0, eta]*A;
        Nzeta = [0, 0, 1, xi]*A;
        
        D = [Neta; Nzeta];
    
        J = D*nodos;

        j(i) = abs(det(J));
        mult = j(i)*w(i);

        Area = Area + mult;
end% i
Area_iso = 4;

dis = j*Area_iso/Area;
dis = min(dis);

end

function [stretch] = alarge(nodos)
    nnod = size(nodos,1);
    L = 0;

    for i = 1:nnod-1
        for j = i+1:nnod
            L(end+1) = norm(nodos(i,:)-nodos(j,:));        
        end
    end
    L(1) = [];

    a = max(L);
    b = min(L);


    stretch = b/a*sqrt(2);
end