function u = solve_FEM(conn, x, y, outb_node, ptop_node, pbottom_node, ...
                             d_elm, e0, er, V1, V2)
    M = size(conn,1); N = length(x);
    eps_e = ones(M,1)*e0; eps_e(d_elm) = er*e0;
    
    A = sparse(N,N); b = sparse(N,1);
    for e = 1:M
        [Ae, ~] = assemble_S(e, x, y, conn, eps_e(e), eps_e(e), 0);
        for i = 1:3
            for j = 1:3
                A(conn(e,i),conn(e,j)) = A(conn(e,i),conn(e,j)) + Ae(i,j);
            end
        end
    end
    
    % Dirichlet: u=0 su bordo esterno, V1/V2 sulle piastre
    BCnodes  = [outb_node; ptop_node; pbottom_node];
    BCvalues = [zeros(size(outb_node)); V1*ones(size(ptop_node)); ...
                V2*ones(size(pbottom_node))];
    nodes = 1:N; nodes(BCnodes) = 0;
    Atemp = speye(N,N); Atemp(find(nodes),:) = 0;
    A(BCnodes,:) = 0; A = A+Atemp; b(BCnodes) = BCvalues;
    
    u = A\b;
    u = full(u);
end
