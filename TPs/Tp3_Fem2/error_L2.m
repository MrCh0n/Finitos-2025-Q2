function [error_L2_abs, error_L2_rel] = error_L2(nodos, elems, P_num, func_analitica)
    % Retorna la norma L2 del error (absoluto) y el error relativo.
    
    nelem = size(elems, 1);
    integral_error2 = 0;
    integral_ana2 = 0;

    % Puntos y pesos de Gauss para Q4 (regla de 2x2)
    puntos = [-1/sqrt(3), 1/sqrt(3)];
    pesos = [1, 1];

    for i = 1:nelem
        nodoid = elems(i, :);
        coord = nodos(nodoid, :);
        P_elem = P_num(nodoid);

        for j = 1:2
            for k = 1:2
                xi = puntos(j);
                eta = puntos(k);
                peso_q = pesos(j) * pesos(k);

                % Funciones de forma N y derivadas locales
                N = 1/4 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
                dN = 1/4 * [-(1-eta),  (1-eta), (1+eta), -(1+eta);
                            -(1-xi), -(1+xi), (1+xi),  (1-xi)];

                % Jacobiano y determinante
                J = dN * coord;
                detJ = det(J);

                % Coordenada global y evaluación escalar
                coord_q = N * coord;
                y_q = coord_q(2);

                P_q_num = N * P_elem;
                P_q_ana = func_analitica(y_q);

                % Acumulación de integrales mediante cuadratura
                integral_error2 = integral_error2 + peso_q * detJ * (P_q_ana - P_q_num)^2;
                integral_ana2 = integral_ana2 + peso_q * detJ * (P_q_ana)^2;
            end
        end
    end

    % Aplicación de raíces para obtener normas
    error_L2_abs = sqrt(integral_error2);
    norma_L2_ana = sqrt(integral_ana2);

    % Cálculo del error relativo con prevención de singularidad
    if norma_L2_ana > 1e-12
        error_L2_rel = error_L2_abs / norma_L2_ana;
    else
        error_L2_rel = 0;
    end
end