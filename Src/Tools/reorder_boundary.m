function B_ordered = reorder_boundary(B_nodes, x, y)
    xB = x(B_nodes);
    yB = y(B_nodes);

    xmin = min(xB); xmax = max(xB);
    ymin = min(yB); ymax = max(yB);
    tol  = 1e-10;

    bot_idx = find(abs(yB - ymin) < tol);
    top_idx = find(abs(yB - ymax) < tol);
    lft_idx = find(abs(xB - xmin) < tol);
    rgt_idx = find(abs(xB - xmax) < tol);

    [~, o] = sort(xB(bot_idx));           bot_idx = bot_idx(o);
    [~, o] = sort(yB(rgt_idx));           rgt_idx = rgt_idx(o);
    [~, o] = sort(xB(top_idx), 'descend'); top_idx = top_idx(o);
    [~, o] = sort(yB(lft_idx), 'descend'); lft_idx = lft_idx(o);

    % Rimuove angoli duplicati
    ordered_idx = [bot_idx; rgt_idx(2:end); top_idx(2:end); lft_idx(2:end-1)];
    B_ordered   = B_nodes(ordered_idx);
end