function T = assemble_T(B_nodes, x, y)
    B_ord = reorder_boundary(B_nodes, x, y);
    n     = length(B_ord);

    % Mappa nodo globale -> indice locale in B_nodes
    global2local = zeros(max(B_nodes), 1);
    for k = 1:length(B_nodes)
        global2local(B_nodes(k)) = k;
    end

    T = zeros(length(B_nodes));

    for k = 1:n
        n1 = B_ord(k);
        n2 = B_ord(mod(k, n) + 1);
        he = sqrt((x(n2)-x(n1))^2 + (y(n2)-y(n1))^2);

        i1 = global2local(n1);
        i2 = global2local(n2);

        T(i1,i1) = T(i1,i1) + he/3;
        T(i1,i2) = T(i1,i2) + he/6;
        T(i2,i1) = T(i2,i1) + he/6;
        T(i2,i2) = T(i2,i2) + he/3;
    end
end
