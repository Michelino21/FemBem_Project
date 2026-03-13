function [] = error_FEM(Lx,Ly,delx,dely,e0,er,V1,V2,u_FEM_BEM,t_FEM_BEM,folder_name,SALVATAGGIO)
    disp("------- ERRORI -------")
    %% === STEP 1: SETUP
    % Calcola C di riferimento

    margin_BEM_coeff = 10;
    margine = margin_BEM_coeff * Ly;     
    Lx_ref = Lx + 2*margine;
    Ly_ref = Ly + 2*margine;

    [conn_r, x_r, y_r, ~, ~, outb_r, ptop_r, pbottom_r, d_elm_r] = ...
        meshgen_cap(Lx_ref, Ly_ref, Lx, Ly, delx, dely);
    [C_ref, ~] = compute_capacitance(conn_r, x_r, y_r, d_elm_r, e0, er, u_FEM_BEM, V1, V2);

    size(u_FEM_BEM)
    
    margini_list = [10,20,30,40,50];
    margini =  margini_list * Ly;
    n_test  = length(margini);
    
    err_L2  = zeros(n_test, 1);
    err_C   = zeros(n_test, 1);
    t_fem   = zeros(n_test, 1);
    N_nodes = zeros(n_test, 1);
    
    %% === STEP 2: FEM PURO CON MARGINI CRESCENTI ===
    
    for k = 1:n_test
        disp(['--- k_fem: ' num2str(k) ' , margin: ' num2str(margini(k)) '---']);
        Lx_k = Lx + 2*margini(k);
        Ly_k = Ly + 2*margini(k);
        
        tic
        % Genera mesh e risolvi FEM puro
        [conn_k, x_k, y_k, ~, ~, outb_k, ptop_k, pbottom_k, d_elm_k] = ...
            meshgen_cap(Lx_k, Ly_k, Lx, Ly, delx, dely);
        
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
        
        % disp(['Margine = ' num2str(margini(k)/Ly) '*Ly  |  ' ...
              %'N = ' num2str(N_nodes(k)) '  |  ' ...
              %'t = ' num2str(t_fem(k)) ' s  |  ' ...
              %'err_L2 = ' num2str(err_L2(k)) '  |  ' ...
              %'err_C = ' num2str(err_C(k)*100) '%'])
    end
    
    %% === STEP 3: PLOT RISULTATI ===
    
    margini_label = margini_list;
    
    % --- Errore L2 vs margine ---
    figure('Color','w');
    semilogy(margini_label, err_L2, 'bo-', 'LineWidth', 2, 'MarkerSize', 8)
    xlabel('Domain margin / Ly');
    ylabel('Relative L^2 error on u');
    title('Pure FEM accuracy vs domain size (reference: FEM-BEM)');
    grid on
    
    salva_fig('L2_norm.png', folder_name, SALVATAGGIO);
    
    % --- Errore capacitanza vs margine ---
    %figure('Color','w');
    %semilogy(margini_label, err_C*100, 'rs-', 'LineWidth', 2, 'MarkerSize', 8)
    %xlabel('Domain margin / Ly');
    %ylabel('Relative error on C (%)');
    %title('Capacitance error vs domain size');
    %grid on
    
    % --- Tempo vs N nodi ---
    figure('Color','w');
    hold on
    plot(margini_list, t_fem, 'bs-', 'LineWidth', 2, 'MarkerSize', 8, ...
         'DisplayName', 'Pure FEM')
    plot(length(x_r), t_FEM_BEM, 'r*', 'MarkerSize', 12, 'LineWidth', 2, ...
         'DisplayName', 'FEM-BEM')
    xlabel('Number of nodes N');
    ylabel('Wall-clock time (s)');
    title('Computational cost comparison');
    legend; grid on

    salva_fig('time.png', folder_name, SALVATAGGIO);
    
    % --- Tabella riassuntiva ---
    fprintf('\n%-15s %-10s %-12s %-10s\n', ...
        'Margin/Ly', 'N nodes', 'err_L2 (%)', 'time (s)'); %'err_C (%)'
    fprintf('%-15s %-10s %-12s %-10s\n', ...
        '---------', '-------', '------', '--------');
    for k = 1:n_test
        fprintf('%-15.0f %-10d %-12.2e %-10.3f\n', ...
            margini_label(k), N_nodes(k), err_L2(k), t_fem(k)); %err_C(k)*100
    end
    fprintf('%-15s %-10d %-12s %-10.3f\n', ...
        'FEM-BEM(ref)', length(x_r), '--', t_FEM_BEM);
end