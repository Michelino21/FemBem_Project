function [] = error_FEM(dx,dy,delx,dely,e0,er,V1,V2,u_FEM_BEM,t_FEM_BEM)
    disp("------- ERRORI -------")
    %% === STEP 1: SETUP
    % Calcola C di riferimento
    margine = 3 * dy;     
    Lx_ref = dx + 2*margine;
    Ly_ref = dy + 2*margine;

    [conn_r, x_r, y_r, ~, ~, outb_r, ptop_r, pbottom_r, d_elm_r] = ...
        meshgen_cap(Lx_ref, Ly_ref, dx, dy, delx, dely);
    [C_ref, ~] = compute_capacitance(conn_r, x_r, y_r, d_elm_r, e0, er, u_FEM_BEM, V1, V2);

    size(u_FEM_BEM)
    
    
    margini = [50,100,150,200] * dy;
    n_test  = length(margini);
    
    err_L2  = zeros(n_test, 1);
    err_C   = zeros(n_test, 1);
    t_fem   = zeros(n_test, 1);
    N_nodes = zeros(n_test, 1);
    
    %% === STEP 2: FEM PURO CON MARGINI CRESCENTI ===
    
    for k = 1:n_test
        disp(['--- k_fem: ' num2str(k) ' , margin: ' num2str(margini(k)) '---']);
        Lx_k = dx + 2*margini(k);
        Ly_k = dy + 2*margini(k);
        
        tic
        % Genera mesh e risolvi FEM puro
        [conn_k, x_k, y_k, ~, ~, outb_k, ptop_k, pbottom_k, d_elm_k] = ...
            meshgen_cap(Lx_k, Ly_k, dx, dy, delx, dely);
        
        u_fem_k = solve_FEM(conn_k, x_k, y_k, outb_k, ptop_k, pbottom_k, ...
                                  d_elm_k, e0, er, V1, V2);
        t_fem(k) = toc;
        N_nodes(k) = length(x_k);
        
        % Errore L2 sui nodi interni comuni (dentro Omega_int del riferimento)
        % Interpola la soluzione FEM sul riferimento
        F_interp = scatteredInterpolant(x_k, y_k, u_fem_k, 'linear', 'nearest');
        u_fem_on_ref = F_interp(x_r, y_r);
        
        % Nodi interni (escludi bordo e piastre)
        inner = setdiff(1:length(x_r), [outb_r; ptop_r; pbottom_r]);

        size(inner)
        
        diff_u = u_fem_on_ref(inner) - u_FEM_BEM(inner);
        err_L2(k) = sqrt(sum(diff_u.^2) / sum(u_FEM_BEM(inner).^2));
        
        % Errore sulla capacitanza
        [C_k, ~] = compute_capacitance(conn_k, x_k, y_k, d_elm_k, e0, er, u_fem_k, V1, V2);
        err_C(k) = abs(C_k - C_ref) / abs(C_ref);
        
        % disp(['Margine = ' num2str(margini(k)/dy) '*dy  |  ' ...
              %'N = ' num2str(N_nodes(k)) '  |  ' ...
              %'t = ' num2str(t_fem(k)) ' s  |  ' ...
              %'err_L2 = ' num2str(err_L2(k)) '  |  ' ...
              %'err_C = ' num2str(err_C(k)*100) '%'])
    end
    
    %% === STEP 3: PLOT RISULTATI ===
    
    margini_label = margini/dy;
    
    % --- Errore L2 vs margine ---
    figure('Color','w');
    semilogy(margini_label, err_L2, 'bo-', 'LineWidth', 2, 'MarkerSize', 8)
    xlabel('Domain margin / dy');
    ylabel('Relative L^2 error on u');
    title('Pure FEM accuracy vs domain size (reference: FEM-BEM)');
    grid on
    
    
    % --- Errore capacitanza vs margine ---
    %figure('Color','w');
    %semilogy(margini_label, err_C*100, 'rs-', 'LineWidth', 2, 'MarkerSize', 8)
    %xlabel('Domain margin / dy');
    %ylabel('Relative error on C (%)');
    %title('Capacitance error vs domain size');
    %grid on
    
    % --- Tempo vs N nodi ---
    figure('Color','w');
    hold on
    plot(N_nodes, t_fem, 'bs-', 'LineWidth', 2, 'MarkerSize', 8, ...
         'DisplayName', 'Pure FEM')
    plot(length(x_r), t_FEM_BEM, 'r*', 'MarkerSize', 12, 'LineWidth', 2, ...
         'DisplayName', 'FEM-BEM')
    xlabel('Number of nodes N');
    ylabel('Wall-clock time (s)');
    title('Computational cost comparison');
    legend; grid on
    
    % --- Tabella riassuntiva ---
    fprintf('\n%-15s %-10s %-12s %-10s\n', ...
        'Margin/dy', 'N nodes', 'err_L2 (%)', 'time (s)'); %'err_C (%)'
    fprintf('%-15s %-10s %-12s %-10s\n', ...
        '---------', '-------', '------', '--------');
    for k = 1:n_test
        fprintf('%-15.0f %-10d %-12.2e %-10.3f\n', ...
            margini_label(k), N_nodes(k), err_L2(k), t_fem(k)); %err_C(k)*100
    end
    fprintf('%-15s %-10d %-12s %-10.3f\n', ...
        'FEM-BEM(ref)', length(x_r), '--', t_FEM_BEM);
end